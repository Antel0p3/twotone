#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
// #include <malloc.h>
// #include <x86intrin.h>
// #include <windows.h>
#include "tool.h"

void rw_args_init(rw_args *irw_args, int64_t line_read_times, int64_t read_cell_num, int64_t line_cell_num, int64_t extra_cell_num, 
    int64_t lineK_num, int64_t *file_read_points, int64_t *file_write_nums, __uint128_t **src_array, __uint128_t **res_array, 
    FILE *src_file, FILE **res_file, sem_t *read_done, sem_t *cal_done, sem_t *write_done)
{
    irw_args->line_read_times = line_read_times;
    irw_args->read_cell_num = read_cell_num;
    irw_args->line_cell_num = line_cell_num;
    irw_args->extra_cell_num = extra_cell_num;
    irw_args->lineK_num = lineK_num;
    irw_args->file_read_points = file_read_points;
    irw_args->file_write_nums = file_write_nums;
    irw_args->src_array = src_array;
    irw_args->res_array = res_array;
    irw_args->src_file = src_file;
    irw_args->res_file = res_file;
    irw_args->read_done = read_done;
    irw_args->cal_done = cal_done;
    irw_args->write_done = write_done;
}

void cal_args_init(cal_args *ical_args, int64_t file_size, int64_t line_read_times, int64_t read_cell_num, int64_t lineK_num, 
    int64_t index, int64_t *file_write_nums, int64_t *res_extras, int64_t *block_offests, __uint128_t **src_array, 
    __uint128_t **res_array, sem_t *read_done, sem_t *cal_done, sem_t *write_done)
{
    ical_args->file_size = file_size;
    ical_args->line_read_times = line_read_times;
    ical_args->read_cell_num = read_cell_num;
    ical_args->lineK_num = lineK_num;
    ical_args->index = index;
    ical_args->file_write_nums = file_write_nums;
    ical_args->res_extras = res_extras;
    ical_args->block_offests = block_offests;
    ical_args->src_array = src_array;
    ical_args->res_array = res_array;
    ical_args->read_done = read_done;
    ical_args->cal_done = cal_done;
    ical_args->write_done = write_done;
}

void time_counter_init(time_counter *time_c)
{
    // QueryPerformanceFrequency(&time_c->freq);
}

void time_counter_begin(time_counter *time_c)
{
    gettimeofday(&(time_c->begin_time), NULL);
    // printf("begin: %ld\n", time_c->begin_time.tv_sec);
    // QueryPerformanceCounter(&time_c->begin_time);
}

void time_counter_end(time_counter *time_c)
{
    gettimeofday(&(time_c->end_time), NULL);
    // printf("end: %ld\n", time_c->end_time.tv_sec);
    // QueryPerformanceCounter(&time_c->end_time);
}

double get_run_time(time_counter time_c)
{
    return (double)(time_c.end_time.tv_sec - time_c.begin_time.tv_sec) * 1000 + (double)(time_c.end_time.tv_usec - time_c.begin_time.tv_usec) / 1000;
    // return MICRO * (double)(time_c.end_time.QuadPart - time_c.begin_time.QuadPart) / (double)time_c.freq.QuadPart;
}

// uint64_t rdtscp() {
//     return __rdtsc();
// }

void *aligned_alloc(size_t _Size, size_t _Alignment)
{
    void * pointer;
	posix_memalign(&pointer, _Alignment, _Size);
	return pointer;
}

void aligned_free(void *_Memory)
{
    free(_Memory);
}

void printbit(void *target, int length)
{
    unsigned char *inner_target = (unsigned char *)target;
    for (int i = 0; i < length; i++)
        printf("%02x ", *(inner_target + i));
    printf("\n");
}

void shiftarray(__uint128_t **arg, long lineK_num, long shift_cell_num, long read_cell_num)
{
    for (int i = 0; i < lineK_num; i++)
        memcpy(arg[i], arg[i] + read_cell_num, shift_cell_num * BYTE_PERCELL);
}

__uint128_t xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift)
{
    __uint128_t arg_res = 0;
    for (int i = 0; i < lineK_num; i++)
        arg_res ^=  arg[i][line_shift + i * row_shift];
    return arg_res;
}

__uint128_t based_xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift, __uint128_t based)
{
    __uint128_t arg_res = based;
    for (int i = 0; i < lineK_num; i++)
        arg_res ^=  arg[i][line_shift + i * row_shift];
    return arg_res;
}

__uint128_t simd_xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift)
{
    __m128i res = _mm_setzero_si128();
    for (int i = 0; i < lineK_num; i++) {
        __m128i tmp = _mm_loadu_si128((__m128i *)&arg[i][line_shift + i *row_shift]);
        res = _mm_xor_si128(res, tmp);
    }
    __uint128_t ret;
    _mm_storeu_si128((__m128i *)&ret, res);
    return ret;
}

__uint128_t simd_based_xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift, __uint128_t based)
{
    __m128i res = _mm_loadu_si128((__m128i *)&based);
    for (int i = 0; i < lineK_num; i++) {
        __m128i tmp = _mm_loadu_si128((__m128i *)&arg[i][line_shift + i *row_shift]);
        res = _mm_xor_si128(res, tmp);
    }
    __uint128_t ret;
    _mm_storeu_si128((__m128i *)&ret, res);
    return ret;
}

void *cal_parallel(void *arg)
{
    cal_args *pcal_args = (cal_args *)arg;
    int64_t file_size = pcal_args->file_size;
    int64_t line_read_times = pcal_args->line_read_times;
    int64_t read_cell_num = pcal_args->read_cell_num;
    int64_t lineK_num = pcal_args->lineK_num;
    int64_t index = pcal_args->index;
    int64_t *file_write_nums = pcal_args->file_write_nums;
    int64_t *res_extras = pcal_args->res_extras;
    int64_t *block_offests = pcal_args->block_offests;
    alignas(CELLSIZE) __uint128_t **src_array = pcal_args->src_array;
    alignas(CELLSIZE) __uint128_t **res_array = pcal_args->res_array;
    sem_t *read_done = pcal_args->read_done;
    sem_t *cal_done = pcal_args->cal_done;
    sem_t *write_done = pcal_args->write_done;
    for (int i = 0; i < line_read_times; i++)
    {
        /* 根据memory中的源数据可以计算出的结果数据量，单位为CELL */
        int64_t need_write_size = read_cell_num < file_write_nums[index] ? read_cell_num : file_write_nums[index];
        /* 根据cache中的源数据可以计算出的结果数据量，单位为CELL */
        int64_t line_max_cache = UPTO_K(MAXCACHE / CELLSIZE, lineK_num) / lineK_num;
        /* 利用完整cache的次数 */
        int64_t cal_times = need_write_size / line_max_cache;

        sem_wait(read_done);

        // 此for循环进行并发xor操作，每次完整利用整个cache
        for(int64_t tt = 0; tt < TEST_TIMES; tt++){
            for (int times = 0; times < cal_times; times++)
                for (int k = 0; k < line_max_cache; k++)
                    res_array[index][times * line_max_cache + k] =
                        simd_xorK(src_array, lineK_num, res_extras[index] + times * line_max_cache + k, block_offests[index]);
            // 此for循环进行并发xor操作，每次部分利用整个cache
            for (int k = 0; k < need_write_size - cal_times * line_max_cache; k++)
                res_array[index][cal_times * line_max_cache + k] =
                    simd_xorK(src_array, lineK_num, res_extras[index] + cal_times * line_max_cache + k, block_offests[index]);
        }

        sem_post(cal_done);
        sem_wait(write_done);
    }
}

// void *rw_parallel(void *arg)
// {
//     rw_args *prw_args = (rw_args *)arg;
//     int64_t line_read_times = prw_args->line_read_times;
//     int64_t read_cell_num = prw_args->read_cell_num;
//     int64_t line_cell_num = prw_args->line_cell_num;
//     int64_t extra_cell_num = prw_args->extra_cell_num;
//     int64_t lineK_num = prw_args->lineK_num;
//     int64_t *file_read_points = prw_args->file_read_points;
//     int64_t *file_write_nums = prw_args->file_write_nums;
//     alignas(CELLSIZE) __uint128_t **src_array = prw_args->src_array;
//     alignas(CELLSIZE) __uint128_t **res_array = prw_args->res_array;
//     FILE *src_file = prw_args->src_file; 
//     FILE **res_file = prw_args->res_file;
//     sem_t *read_done = prw_args->read_done;
//     sem_t *cal_done = prw_args->cal_done;
//     sem_t *write_done = prw_args->write_done;
//     int64_t counter = 0;
//     for (int i = 0; i < line_read_times; i++)
//     {
//         // 每次读取的数据量不能超过read_cell_num或剩余数据line_cell_num - counter
//         int64_t need_read_size = read_cell_num < line_cell_num - counter ? read_cell_num : line_cell_num - counter;
//         // 此for循环进行数据读取操作
//         for (int j = 0; j < lineK_num; j++)
//         {
//             fseek(src_file, file_read_points[j], SEEK_SET);
//             fread(src_array[j] + extra_cell_num, sizeof(__uint128_t), need_read_size, src_file);
//             file_read_points[j] += need_read_size * BYTE_PERCELL;
//         }
//         counter += need_read_size;

//         sem_post_multiple(read_done, (int)lineK_num);
//         for (int j = 0; j < lineK_num; j++)
//             sem_wait(cal_done);

//         /* 根据内存中的源数据可以计算出的结果数据量，单位为CELL */
//         int64_t need_write_size[lineK_num];
//         // 此for循环本轮每个冗余块最多可以计算的数据量，j是第j个冗余块
//         for (int j = 0; j < lineK_num; j++)
//         {
//             need_write_size[j] = read_cell_num < file_write_nums[j] ? read_cell_num : file_write_nums[j];
//             file_write_nums[j] = 0 > file_write_nums[j] - read_cell_num ? 0 : file_write_nums[j] - read_cell_num;
//         }
//         for (int j = 0; j < lineK_num; j++)
//             fwrite(res_array[j], sizeof(__uint128_t), need_write_size[j], res_file[j]);
//         for (int j = 0; j < lineK_num; j++)
//             fflush(res_file[j]);
//         // 此处进行extra_cell的更新
//         shiftarray(src_array, lineK_num, extra_cell_num, read_cell_num);
//         for (int j = 0; j < lineK_num; j++)
//             memset(src_array[j] + extra_cell_num, 0, read_cell_num * BYTE_PERCELL);
//         for (int j = 0; j < lineK_num; j++)
//             memset(res_array[j], 0, read_cell_num * BYTE_PERCELL);
        
//         sem_post_multiple(write_done, (int)lineK_num);
//     }
// }

char *realname(char *name)
{
    int counter = 0;
    char *oldname = name;
    while (*(name++) != '\0')
    {
        if (name[0] == '.' && name[1] == 'b' && name[2] == 'a' && name[3] == 'c' &&
             name[4] == 'k' && name[5] == 'u' && name[6] == 'p' && name[7] == '.' && name[8] == 'v')
            break;
        counter++;
    }
    char *realname = calloc(counter + 2, sizeof(char));
    memcpy(realname, oldname, counter + 1);
    return realname;
}

long getK(char *name)
{
    int counter_start = 0, counter_over = 0;
    int flag = 1;
    while (*(name++) != '\0')
    {
        if (name[0] == '.' && name[1] == 'b' && name[2] == 'a' && name[3] == 'c' &&
             name[4] == 'k' && name[5] == 'u' && name[6] == 'p' && name[7] == '.' && name[8] == 'v')
            break;
        counter_start++;
    }
    while (*(name++) != '\0')
    {
        if (name[0] == '(' && name[1] == 'k' && name[2] == '=')
            break;
        counter_start++;
    }
    counter_over = counter_start;
    char *kname = name + 3;
    while (*(name++) != '\0')
    {
        if (name[0] == ')')
        {
            flag = 0;
            break;
        }
        counter_over++;
    }
    if (flag)
    {
        printf("非backup文件，格式错误！");
        exit(-1);
    }
    char *realk = calloc(counter_over - counter_start - 1, sizeof(char));
    memcpy(realk, kname, counter_over - counter_start - 2);
    return atoi(realk);
}

void encodek3(long line_cell_num, FILE *src_file, FILE *res_file)
{
    __uint128_t *arrayk3[3];
    __uint128_t *res_array;
    arrayk3[0] = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    arrayk3[1] = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    arrayk3[2] = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk3[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk3[1] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk3[2] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 2; i++)
        res_array[i] = arrayk3[0][i] ^ arrayk3[1][i] ^ arrayk3[2][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 2, res_file);
    free(arrayk3[0]);
    free(arrayk3[1]);
    free(arrayk3[2]);

    arrayk3[0] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk3[1] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk3[2] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk3[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk3[1], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk3[2], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num; i++)
        res_array[i] = arrayk3[0][i] ^ arrayk3[1][i] ^ arrayk3[2][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num, res_file);
    free(arrayk3[0]);
    free(arrayk3[1]);
    free(arrayk3[2]);
    free(res_array);

    arrayk3[0] = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    arrayk3[1] = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    arrayk3[2] = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 2, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk3[0] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk3[1] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk3[2], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 2; i++)
        res_array[i] = arrayk3[0][i] ^ arrayk3[1][i] ^ arrayk3[2][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 2, res_file);
    free(arrayk3[0]);
    free(arrayk3[1]);
    free(arrayk3[2]);
    free(res_array);
}

void encodek4(long line_cell_num, FILE *src_file, FILE *res_file)
{
    __uint128_t *arrayk4[4];
    __uint128_t *res_array;
    arrayk4[0] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    arrayk4[1] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    arrayk4[2] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    arrayk4[3] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk4[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[1] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[2] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[3] + 3, sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 3; i++)
        res_array[i] = arrayk4[0][i] ^ arrayk4[1][i] ^ arrayk4[2][i] ^ arrayk4[3][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 3, res_file);
    free(arrayk4[0]);
    free(arrayk4[1]);
    free(arrayk4[2]);
    free(arrayk4[3]);
    free(res_array);

    arrayk4[0] = (__uint128_t *)calloc(line_cell_num * sizeof(__uint128_t), CELLSIZE);
    arrayk4[1] = (__uint128_t *)calloc(line_cell_num * sizeof(__uint128_t), CELLSIZE);
    arrayk4[2] = (__uint128_t *)calloc(line_cell_num * sizeof(__uint128_t), CELLSIZE);
    arrayk4[3] = (__uint128_t *)calloc(line_cell_num * sizeof(__uint128_t), CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num * sizeof(__uint128_t), CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk4[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[1], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[2], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[3], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num; i++)
        res_array[i] = arrayk4[0][i] ^ arrayk4[1][i] ^ arrayk4[2][i] ^ arrayk4[3][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num, res_file);
    free(arrayk4[0]);
    free(arrayk4[1]);
    free(arrayk4[2]);
    free(arrayk4[3]);
    free(res_array);

    arrayk4[0] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    arrayk4[1] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    arrayk4[2] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    arrayk4[3] = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 3, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk4[0] + 3, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[1] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[2] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[3], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 3; i++)
        res_array[i] = arrayk4[0][i] ^ arrayk4[1][i] ^ arrayk4[2][i] ^ arrayk4[3][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 3, res_file);
    free(arrayk4[0]);
    free(arrayk4[1]);
    free(arrayk4[2]);
    free(arrayk4[3]);
    free(res_array);

    arrayk4[0] = (__uint128_t *)calloc(line_cell_num + 6, CELLSIZE);
    arrayk4[1] = (__uint128_t *)calloc(line_cell_num + 6, CELLSIZE);
    arrayk4[2] = (__uint128_t *)calloc(line_cell_num + 6, CELLSIZE);
    arrayk4[3] = (__uint128_t *)calloc(line_cell_num + 6, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 6, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk4[0] + 6, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[1] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[2] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk4[3], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 6; i++)
        res_array[i] = arrayk4[0][i] ^ arrayk4[1][i] ^ arrayk4[2][i] ^ arrayk4[3][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 6, res_file);
    free(arrayk4[0]);
    free(arrayk4[1]);
    free(arrayk4[2]);
    free(arrayk4[3]);
    free(res_array);
}

void encodek5(long line_cell_num, FILE *src_file, FILE *res_file)
{
    __uint128_t *arrayk5[5], *res_array;
    arrayk5[0] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[1] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[2] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[3] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[4] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk5[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[1] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[2] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[3] + 6, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[4] + 8, sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 8; i++)
        res_array[i] = arrayk5[0][i] ^ arrayk5[1][i] ^ arrayk5[2][i] ^ arrayk5[3][i] ^ arrayk5[4][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 8, res_file);
    free(arrayk5[0]);
    free(arrayk5[1]);
    free(arrayk5[2]);
    free(arrayk5[3]);
    free(arrayk5[4]);
    free(res_array);

    arrayk5[0] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[1] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[2] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[3] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[4] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk5[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[1] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[2] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[3] + 3, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[4] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 4; i++)
        res_array[i] = arrayk5[0][i] ^ arrayk5[1][i] ^ arrayk5[2][i] ^ arrayk5[3][i] ^ arrayk5[4][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 4, res_file);
    free(arrayk5[0]);
    free(arrayk5[1]);
    free(arrayk5[2]);
    free(arrayk5[3]);
    free(arrayk5[4]);
    free(res_array);

    arrayk5[0] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk5[1] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk5[2] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk5[3] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk5[4] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk5[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[1], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[2], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[3], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[4], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num; i++)
        res_array[i] = arrayk5[0][i] ^ arrayk5[1][i] ^ arrayk5[2][i] ^ arrayk5[3][i] ^ arrayk5[4][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num, res_file);
    free(arrayk5[0]);
    free(arrayk5[1]);
    free(arrayk5[2]);
    free(arrayk5[3]);
    free(arrayk5[4]);
    free(res_array);

    arrayk5[0] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[1] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[2] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[3] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    arrayk5[4] = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 4, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk5[0] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[1] + 3, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[2] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[3] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[4], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 4; i++)
        res_array[i] = arrayk5[0][i] ^ arrayk5[1][i] ^ arrayk5[2][i] ^ arrayk5[3][i] ^ arrayk5[4][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 4, res_file);
    free(arrayk5[0]);
    free(arrayk5[1]);
    free(arrayk5[2]);
    free(arrayk5[3]);
    free(arrayk5[4]);
    free(res_array);

    arrayk5[0] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[1] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[2] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[3] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    arrayk5[4] = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 8, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk5[0] + 8, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[1] + 6, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[2] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[3] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk5[4], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 8; i++)
        res_array[i] = arrayk5[0][i] ^ arrayk5[1][i] ^ arrayk5[2][i] ^ arrayk5[3][i] ^ arrayk5[4][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 8, res_file);
    free(arrayk5[0]);
    free(arrayk5[1]);
    free(arrayk5[2]);
    free(arrayk5[3]);
    free(arrayk5[4]);
    free(res_array);
}

void encodek6(long line_cell_num, FILE *src_file, FILE *res_file)
{
    __uint128_t *arrayk6[6], *res_array;
    arrayk6[0] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[1] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[2] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[3] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[4] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[5] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk6[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[1] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[2] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[3] + 6, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[4] + 8, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[5] + 10, sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 10; i++)
        res_array[i] = arrayk6[0][i] ^ arrayk6[1][i] ^ arrayk6[2][i] ^ arrayk6[3][i] ^ arrayk6[4][i] ^ arrayk6[5][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 10, res_file);
    free(arrayk6[0]);
    free(arrayk6[1]);
    free(arrayk6[2]);
    free(arrayk6[3]);
    free(arrayk6[4]);
    free(arrayk6[5]);
    free(res_array);

    arrayk6[0] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[1] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[2] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[3] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[4] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[5] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk6[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[1] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[2] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[3] + 3, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[4] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[5] + 5, sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 5; i++)
        res_array[i] = arrayk6[0][i] ^ arrayk6[1][i] ^ arrayk6[2][i] ^ arrayk6[3][i] ^ arrayk6[4][i] ^ arrayk6[5][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 5, res_file);
    free(arrayk6[0]);
    free(arrayk6[1]);
    free(arrayk6[2]);
    free(arrayk6[3]);
    free(arrayk6[4]);
    free(arrayk6[5]);
    free(res_array);

    arrayk6[0] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk6[1] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk6[2] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk6[3] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk6[4] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    arrayk6[5] = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk6[0], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[1], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[2], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[3], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[4], sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[5], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num; i++)
        res_array[i] = arrayk6[0][i] ^ arrayk6[1][i] ^ arrayk6[2][i] ^ arrayk6[3][i] ^ arrayk6[4][i] ^ arrayk6[5][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num, res_file);
    free(arrayk6[0]);
    free(arrayk6[1]);
    free(arrayk6[2]);
    free(arrayk6[3]);
    free(arrayk6[4]);
    free(arrayk6[5]);
    free(res_array);

    arrayk6[0] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[1] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[2] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[3] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[4] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    arrayk6[5] = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 5, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk6[0] + 5, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[1] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[2] + 3, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[3] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[4] + 1, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[5], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 5; i++)
        res_array[i] = arrayk6[0][i] ^ arrayk6[1][i] ^ arrayk6[2][i] ^ arrayk6[3][i] ^ arrayk6[4][i] ^ arrayk6[5][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 5, res_file);
    free(arrayk6[0]);
    free(arrayk6[1]);
    free(arrayk6[2]);
    free(arrayk6[3]);
    free(arrayk6[4]);
    free(arrayk6[5]);
    free(res_array);

    arrayk6[0] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[1] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[2] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[3] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[4] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    arrayk6[5] = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 10, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk6[0] + 10, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[1] + 8, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[2] + 6, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[3] + 4, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[4] + 2, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[5], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 10; i++)
        res_array[i] = arrayk6[0][i] ^ arrayk6[1][i] ^ arrayk6[2][i] ^ arrayk6[3][i] ^ arrayk6[4][i] ^ arrayk6[5][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 10, res_file);
    free(arrayk6[0]);
    free(arrayk6[1]);
    free(arrayk6[2]);
    free(arrayk6[3]);
    free(arrayk6[4]);
    free(arrayk6[5]);
    free(res_array);

    arrayk6[0] = (__uint128_t *)calloc(line_cell_num + 15, CELLSIZE);
    arrayk6[1] = (__uint128_t *)calloc(line_cell_num + 15, CELLSIZE);
    arrayk6[2] = (__uint128_t *)calloc(line_cell_num + 15, CELLSIZE);
    arrayk6[3] = (__uint128_t *)calloc(line_cell_num + 15, CELLSIZE);
    arrayk6[4] = (__uint128_t *)calloc(line_cell_num + 15, CELLSIZE);
    arrayk6[5] = (__uint128_t *)calloc(line_cell_num + 15, CELLSIZE);
    res_array = (__uint128_t *)calloc(line_cell_num + 15, CELLSIZE);
    fseek(src_file, 0, SEEK_SET);
    fread(arrayk6[0] + 15, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[1] + 12, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[2] + 9, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[3] + 6, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[4] + 3, sizeof(__uint128_t), line_cell_num, src_file);
    fread(arrayk6[5], sizeof(__uint128_t), line_cell_num, src_file);
    fseek(src_file, 0, SEEK_SET);
    for (int i = 0; i < line_cell_num + 15; i++)
        res_array[i] = arrayk6[0][i] ^ arrayk6[1][i] ^ arrayk6[2][i] ^ arrayk6[3][i] ^ arrayk6[4][i] ^ arrayk6[5][i];
    fwrite(res_array, sizeof(__uint128_t), line_cell_num + 15, res_file);
    free(arrayk6[0]);
    free(arrayk6[1]);
    free(arrayk6[2]);
    free(arrayk6[3]);
    free(arrayk6[4]);
    free(arrayk6[5]);
    free(res_array);
}

void espf(char *fna, int64_t ver, int64_t fsz, int64_t k, double ct, int64_t sum, int64_t test_times)
{
    FILE *wf = fopen(fna, "a");
    char temp[2048];
    if (wf == NULL)
    {
        printf("espf fopen fail\n");
        exit(-1);
    }
    // encode version, file size, arg k, cost time
    int64_t file_size = fsz / BYTESIZE;
    int64_t throughout = (fsz + sum * CELLSIZE) / BYTESIZE * test_times;
    sprintf(temp, "encode_ver %lld %lld %lld %.2lf %lld\n", ver, file_size, k, ct, throughout);
    fwrite(temp, sizeof(char), strlen(temp), wf);
    fclose(wf);
}

void dspf(char *fna, int64_t ver, int64_t fsz, int64_t k, double ct, int64_t sum, int64_t test_times)
{
    FILE *wf = fopen(fna, "a");
    char temp[2048];
    if (wf == NULL)
    {
        printf("dspf fopen fail\n");
        exit(-1);
    }
    // decode version, file size, arg k, cost time
    int64_t file_size = fsz / BYTESIZE;
    int64_t throughout = (sum * CELLSIZE) / BYTESIZE * test_times;
    sprintf(temp, "decode_ver %lld %lld %lld %.0lf %lld\n", ver, file_size, k, ct, throughout);
    fwrite(temp, sizeof(char), strlen(temp), wf);
    fclose(wf);
}