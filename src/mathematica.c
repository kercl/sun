#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "WolframLibrary.h"

#include "keylist.h"
#include "irrep.h"
#include "int_gt.h"

struct key_list_node *mngd_basis_list = NULL;

DLLEXPORT void manage_SUNGTBasis(WolframLibraryData libData, mbool mode, mint id)
{
    if(mode == 0) {
      printf("Creating basis\n");
      struct gt_tree *new_tree = malloc(sizeof(struct gt_tree));
      new_tree->num_patterns = 0;
      mngd_basis_list = kl_push(mngd_basis_list, id, (void*)new_tree);
    } else {
      printf("Releasing basis\n");

      struct gt_tree *tree;
      mngd_basis_list = kl_unlink(mngd_basis_list, (int)id, (void**)&tree);

      if(tree == NULL)
        return;

      if(tree->num_patterns > 0)
        gt_free_tree(tree, 1);
      free(tree);
    }
}

DLLEXPORT mint WolframLibrary_getVersion() {
  return WolframLibraryVersion;
}

DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
  int err = (*libData->registerLibraryExpressionManager)("SUNGTBasis", &manage_SUNGTBasis);
  
  return err;
}

DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) {
  (*libData->unregisterLibraryExpressionManager)("SUNGTBasis");
}

void get_top_row_from_args(WolframLibraryData libData, MArgument *Args, gt_int_t **top_row, size_t *len) {
  MTensor arg = MArgument_getMTensor(Args[0]);

  *len = libData->MTensor_getFlattenedLength(arg);
  mint* data = libData->MTensor_getIntegerData(arg);

  *top_row = malloc(sizeof(gt_int_t) * (*len));
  for(int i = 0; i < *len; ++i)
    (*top_row)[i] = (gt_int_t)data[i];
}

DLLEXPORT int ml_dimension_irrep(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
  gt_int_t *top_row;
  size_t len;

  get_top_row_from_args(libData, Args, &top_row, &len);

  MArgument_setInteger(Res, (int)gt_num_of_patterns(top_row, len));

  free(top_row);
  return LIBRARY_NO_ERROR;
}

DLLEXPORT int ml_init_basis(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
  struct gt_tree *basis;

  int basis_id = MArgument_getInteger(Args[1]);
  int ret = kl_find(mngd_basis_list, basis_id, (void**)&basis);

  if(ret == 0) {
    libData->Message("No instance of this basis."); 
    return LIBRARY_FUNCTION_ERROR; 
  }
  
  if(basis->num_patterns > 0) {
    libData->Message("Basis already initialized."); 
    return LIBRARY_FUNCTION_ERROR; 
  }

  gt_int_t *top_row;
  gt_int_t *patterns;
  size_t num_patterns, len;
    
  get_top_row_from_args(libData, Args, &top_row, &len);
  gt_generate_all(&patterns, &num_patterns, top_row, len);
  gt_list_to_tree(basis, patterns, num_patterns, len);
  free(top_row);

  return LIBRARY_NO_ERROR;
}

DLLEXPORT int ml_cartan_matrix(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
  struct gt_tree *basis;

  int basis_id = MArgument_getInteger(Args[0]),
      l = MArgument_getInteger(Args[1]);

  if(kl_find(mngd_basis_list, basis_id, (void**)&basis) == 0) {
    libData->Message("No instance of this basis."); 
    return LIBRARY_FUNCTION_ERROR; 
  }

  MTensor out;
  mint out_type = MType_Integer;
  mint out_rank = 1;
  mint out_dims[] = {basis->num_patterns};
  mint* out_data;
  int err;

  err = libData->MTensor_new( out_type, out_rank, out_dims, &out );
  out_data = libData->MTensor_getIntegerData(out);

  mat_int_t *diagonal = malloc(sizeof(mat_int_t) * basis->num_patterns);
  csa_generator_diag_from_gt(basis, l, diagonal);

  for(int i = 0; i < basis->num_patterns; i++)
    out_data[i] = diagonal[i];

  free(diagonal);
  MArgument_setMTensor(Res, out);
  return LIBRARY_NO_ERROR;
}

DLLEXPORT int ml_root_lowering(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
  struct gt_tree *basis;

  int basis_id = MArgument_getInteger(Args[0]),
      l = MArgument_getInteger(Args[1]);

  if(kl_find(mngd_basis_list, basis_id, (void**)&basis) == 0) {
    libData->Message("No instance of this basis."); 
    return LIBRARY_FUNCTION_ERROR; 
  }

  mat_int_t *numerators, *denominators;
  size_t *rows, *cols, n_entries;
  n_entries = lowering_operator_from_gt(basis, l, &numerators, &denominators, &rows, &cols);

  MTensor out;
  mint out_type = MType_Integer;
  mint out_rank = 2;
  mint out_dims[] = {4, n_entries};
  mint* out_data;
  int err;

  err = libData->MTensor_new( out_type, out_rank, out_dims, &out );
  out_data = libData->MTensor_getIntegerData(out);

  for(int i = 0; i < n_entries; i++) {
    out_data[i] = numerators[i];
    out_data[n_entries + i] = denominators[i];
    out_data[2 * n_entries + i] = rows[i] + 1;
    out_data[3 * n_entries + i] = cols[i] + 1;
  }

  free(numerators);
  free(denominators);
  free(rows);
  free(cols);

  MArgument_setMTensor(Res, out);
  return LIBRARY_NO_ERROR;
}
