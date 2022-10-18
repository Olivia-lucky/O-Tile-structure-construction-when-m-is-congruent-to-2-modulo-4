#ifndef PERMUANDCOMBI2_H
#define PERMUANDCOMBI2_H


#include <cstdio>
#include <cstdlib>
#include <cstring>


class permuandcombi2 {
public:
    //private:
    int st;  // starting point 0 or 1
    int ed;  // st, ..., ed-1
    int g;  // 0 < g <= ed-s 
    int* len;
    int* arr;

public:
    permuandcombi2(int s_st, int s_ed, int s_g) {
        st = s_st;
        ed = s_ed;
        g = s_g;
        len = (int*)malloc(2 * g * sizeof(int));
        arr = len + g;
        for (int i = 0; i < g; i++) {
            len[i] = st + i;
        }
        memcpy(arr, len, g * sizeof(int));
    }
    ~permuandcombi2() {
        free(len);
    }

    int next() {
        if (!next_arr()) {
            if (!next_len()) {
                return 0;
            }
            memcpy(arr, len, g * sizeof(int));
            return 1;
        }
        return 1;
    }
    //private:
    int next_len() {
        int i = g - 1;
        while ((i >= 0) && len[i] == ed - g + i) {
            i--;
        }
        if (i < 0) return 0;
        len[i]++;
        for (int j = i + 1; j < g; j++) {
            len[j] = len[i] + j - i;
        }
        return 1;
    }

    int next_arr() {
        int j = g - 2;
        while ((j >= 0) && arr[j] > arr[j + 1]) {
            j--;
        }
        if (j < 0) return 0;
        int k = g - 1;
        while (arr[j] > arr[k]) {
            k--;
        }

        // swap arr[j and k]
        int temp = arr[k];
        arr[k] = arr[j];
        arr[j] = temp;

        int start = j + 1;
        int end = g - 1;
        while (start < end) {
            temp = arr[start];
            arr[start] = arr[end];
            arr[end] = temp;

            start++;
            end--;
        }
        return 1;
    }
};



#endif // #ifndef PERMUANDCOMBI2_H



#pragma once


