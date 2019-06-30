/* Copyright (c) 2018 Clemens Kerschbaumer
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdlib.h>

#include "keylist.h"

struct key_list_node* kl_push(struct key_list_node *list, int key, void *data) {
    struct key_list_node *new_node = malloc(sizeof(struct key_list_node));

    new_node->key = key;
    new_node->data = data;
    new_node->next = list;

    return new_node;
}

int kl_find(struct key_list_node *list, int key, void **data) {
    while(list) {
        if(list->key == key) {
            *data = list->data;
            return KL_OK;
        }

        list = list->next;
    }
    return KL_MISSING;
}

struct key_list_node* kl_unlink(struct key_list_node *list, int key, void **data) {
    if(list == NULL) {
        *data = NULL;
        return NULL;
    }

    if(list->key == key) {
        *data = list->data;
        struct key_list_node *next = list->next;
        free(list);
        return next;
    }

    struct key_list_node *prev = list, *head = list;
    list = list->next;
    while(list != NULL) {
        if(list->key == key) {
            prev->next = list->next;
            *data = list->data;
            free(list);
            return head;
        }
        prev = prev->next;
        list = list->next;
    }

    *data = NULL;
    return head;
}
