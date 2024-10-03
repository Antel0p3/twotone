#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "tool.c"

// src_file源文件，res_file结果文件(多个写入点，并发写入)
FILE *src_file, **res_file;
// src_array源数组，res_array结果数组
alignas(CELLSIZE) __uint128_t **src_array, **res_array;

int main(int argc, char **argv)
{
    if (argc != 2 && argc != 3)
    {
        printf("please input the ness arg1: filename with the sel arg2: k\n");
        exit(-1);
    }
    // 编码矩阵的列数，也是文件划分后的行数
    int64_t lineK_num = (argc == 2) ? 3 : atoi(argv[2]);
    // 读取源文件src_file，并初始化结果文件res_file
    char *read_file_name = argv[1];
    char *write_file_name = calloc(strlen(read_file_name) + 30, sizeof(char));
    sprintf(write_file_name, "en3_k%ld_%s", lineK_num, read_file_name);
    src_file = fopen(read_file_name, "rb");
    if (!src_file)
    {
        printf("no such file, please input right filename\n");
        exit(-1);
    }
    res_file = calloc(lineK_num, sizeof(FILE *));
    for (int i = 0; i < lineK_num; i++)
        res_file[i] = fopen(write_file_name, "wb");
    // 读取源文件大小
    fseek(src_file, 0, SEEK_END);
    // 文件大小，单位为bit
    int64_t file_size = ftell(src_file) * BYTESIZE;
    // 文件可划分的总CELL个数(已对齐)，先将文件大小对齐到CELLSIZE，再将对应的总CELL个数对齐到K
    int64_t all_cell_num = UPTO_K(UPTO_K(file_size, CELLSIZE) / CELLSIZE, lineK_num);
    // 文件划分后每一行的CELL个数，包含了补零行
    int64_t line_cell_num = all_cell_num / lineK_num;
    // src_array每一行对应的CELL个数(已对齐)
    int64_t read_cell_num = 0;
    // 用于文件过大时，缓存本次计算中所需的上次计算出的部分数据
    int64_t extra_cell_num = (lineK_num - 1) * (lineK_num / 2);
    /* 源文件的读取划分点，用于将src_array的每一行映射到文件对应位置，单位为byte */
    int64_t file_read_points[lineK_num];
    /* 用于记录结果文件的每个文件写入点应该写入的数据量，单位为CELL */ 
    int64_t file_write_nums[lineK_num];
    /* 记录每一个冗余块的偏移单位 */ 
    int64_t block_offests[lineK_num];
    /* 记录每一个冗余块的偏移总和 */
    int64_t block_extras[lineK_num];
    /* xorK异或起始点的绝对偏移量 */
    int64_t res_extras[lineK_num];
    // 用于当文件大小不为CELLSIZE*line_cell_num的整数倍时，避免读取超过改行最大限制的数据
    int64_t counter = 0;
    // 表示中间的冗余块的位置
    int64_t the_middle = (lineK_num - 1) / 2;


    // 文件过大的情况，Memory装不下整个文件，需要多次读取文件内容
    // 将MAXMEMORY对齐到CELLSIZE，再将对应的CELL个数对齐到K，再计算src_array每行的CELL个数
    if (all_cell_num * CELLSIZE > MAXMEMORY)
        read_cell_num = UPTO_K(UPTO_K(MAXMEMORY, CELLSIZE) / CELLSIZE, lineK_num) / lineK_num;
    // 文件不超限的情况，Memory装的下整个文件，无需多次读取文件内容
    // src_array每行的CELL个数就等于文件划分后每一行的CELL个数
    else
        read_cell_num = line_cell_num;
    // 在考虑shift后，读取完文件划分出的一行需要的读取次数
    int64_t line_read_times = UPTO_K(line_cell_num + extra_cell_num, read_cell_num) / read_cell_num;


    // src_array源数组，res_array结果数组的内存初始化
    // src_array的结构为(lineK_num行，read_cell_num + extra_cell_num列)的二维数组
    src_array = (__uint128_t **)aligned_alloc(lineK_num * sizeof(__uint128_t *), CELLSIZE);
    for (int i = 0; i < lineK_num; i++)
        src_array[i] = (__uint128_t *)aligned_alloc((read_cell_num + extra_cell_num) * sizeof(__uint128_t), CELLSIZE);
    // res_array的结构为(lineK_num行，read_cell_num + extra_cell_num列)的二维数组
    res_array = (__uint128_t **)aligned_alloc(lineK_num * sizeof(__uint128_t *), CELLSIZE);
    for (int i = 0; i < lineK_num; i++)
        res_array[i] = (__uint128_t *)aligned_alloc(read_cell_num * sizeof(__uint128_t), CELLSIZE);
    for (int i = 0; i < lineK_num; i++)
        memset(src_array[i], 0, (read_cell_num + extra_cell_num) * BYTE_PERCELL);
    for (int i = 0; i < lineK_num; i++)
        memset(res_array[i], 0, read_cell_num * BYTE_PERCELL);


    // 通过空间局部性加速计算
    // 计算每个冗余块的偏移单位
    for (int i = 0; i < lineK_num; i++)
        block_offests[i] = i - the_middle;
    // 计算每个冗余块的偏移总和
    for (int i = 0; i < lineK_num; i++)
        block_extras[i] = (lineK_num - 1) * labs(block_offests[i]);
    // 计算xorK异或起始点的绝对偏移量
    for (int i = 0; i < lineK_num; i++)
        res_extras[i] = i <= the_middle ? extra_cell_num : extra_cell_num - block_extras[i];
        // 源文件的读取划分点初始化，对于不超限的文件也需要设置file_read_points
    for (int i = 0; i < lineK_num; i++)
        file_read_points[i] = line_cell_num * i * BYTE_PERCELL;
    // 确定每个结果文件的文件写入点应该写入的数据量
    for (int i = 0; i < lineK_num; i++)
        file_write_nums[i] = line_cell_num + block_extras[i];
    // 结果文件的写入划分点初始化，用于并发写入时确定每个写入点位置
    // 每个写入点写入数量不同，但关于零偏移组对称
    int64_t sum = 0;
    for (int64_t i = 0; i < lineK_num; i++)
    {
        fseek(res_file[i], sum * BYTE_PERCELL, SEEK_SET);
        sum += file_write_nums[i];
    }

    time_counter time_c;
    time_counter_init(&time_c);
    time_counter_begin(&time_c);
    // 此for循环进行结果的计算，i是读取的轮次
    for (int i = 0; i < line_read_times; i++)
    {
        // 每次读取的数据量不能超过read_cell_num或剩余数据line_cell_num - counter
        int64_t need_read_size = read_cell_num < line_cell_num - counter ? read_cell_num : line_cell_num - counter;
        // 此for循环进行数据读取操作
        for (int j = 0;  j < lineK_num; j++)
        {
            fseek(src_file, file_read_points[j], SEEK_SET);
            fread(src_array[j] + extra_cell_num, sizeof(__uint128_t), need_read_size, src_file);
            file_read_points[j] += need_read_size * BYTE_PERCELL;
        }
        counter += need_read_size;

        /* 根据内存中的源数据可以计算出的结果数据量，单位为CELL */
        int64_t need_write_size[lineK_num];
        // 此for循环本轮每个冗余块最多可以计算的数据量，j是第j个冗余块
        for (int j = 0; j < lineK_num; j++)
        {
            need_write_size[j] = read_cell_num < file_write_nums[j] ? read_cell_num : file_write_nums[j];
            file_write_nums[j] = 0 > file_write_nums[j] - read_cell_num ? 0 : file_write_nums[j] - read_cell_num;
        }

        /* 根据cache中的源数据可以计算出的结果数据量，单位为CELL */
        int64_t line_max_cache = UPTO_K(MAXCACHE / CELLSIZE, lineK_num) / lineK_num;
        int64_t cal_times = need_write_size[lineK_num - 1] / line_max_cache;
        // 此for循环进行并发xor操作
        for(int64_t tt = 0; tt < TEST_TIMES; tt++){
            for (int times = 0; times < cal_times; times++)
                for (int j = 0; j < lineK_num; j++)
                    for (int k = 0; k < line_max_cache; k++)
                        res_array[j][times * line_max_cache + k] = 
                            simd_xorK(src_array, lineK_num, res_extras[j] + times * line_max_cache + k, block_offests[j]);
            for (int j = 0; j < lineK_num; j++)
                for (int k = 0; k < need_write_size[j] - cal_times * line_max_cache; k++)
                    res_array[j][cal_times * line_max_cache + k] = 
                        simd_xorK(src_array, lineK_num, res_extras[j] + cal_times * line_max_cache + k, block_offests[j]);
        }
        for (int j = 0; j < lineK_num; j++)
            fwrite(res_array[j], sizeof(__uint128_t), need_write_size[j], res_file[j]);

        // 此处进行extra_cell的更新
        shiftarray(src_array, lineK_num, extra_cell_num, read_cell_num);
        for (int j = 0; j < lineK_num; j++)
            memset(src_array[j] + extra_cell_num, 0, read_cell_num * BYTE_PERCELL);
        for (int j = 0; j < lineK_num; j++)
            memset(res_array[j], 0, read_cell_num * BYTE_PERCELL);
    }
    time_counter_end(&time_c);
    // printf("encode ver3, filename = %s, k = %ld, cost time = %lu\n", read_file_name, lineK_num, finish_time-start_time);
    espf("./ente.txt", 3, file_size, lineK_num, get_run_time(time_c), sum, TEST_TIMES);

    for (int i = 0; i < lineK_num; i++)
        aligned_free(src_array[i]);
    aligned_free(src_array);
    for (int i = 0; i < lineK_num; i++)
        aligned_free(res_array[i]);
    aligned_free(res_array);
    fclose(src_file);
    for (int i = 0; i < lineK_num; i++)
        fclose(res_file[i]);
}