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
        printf("please input the filename\n");
        exit(-1);
    }
    // 读取源文件src_file，并初始化结果文件res_file
    char *read_file_name = argv[1];
    int64_t lineK_num = (argc == 2) ? 3 : atoi(argv[2]);
    char *write_file_name = calloc(strlen(read_file_name) + 50, sizeof(char));
    sprintf(write_file_name, "%s_de4", realname(read_file_name));
    src_file = fopen(read_file_name, "rb");
    if (!src_file)
    {
        printf("no such file, please input right filename\n");
        exit(-1);
    }
    // 读取源文件大小
    fseek(src_file, 0, SEEK_END);
    // 文件大小，单位为bit
    int64_t file_size = ftell(src_file) * BYTESIZE;
    // 冗余块的文件大小，单位为CELL
    int64_t all_cell_num = file_size / CELLSIZE;
    // 原始块的文件大小，单位为CELL
    int64_t true_cell_num = all_cell_num;
    res_file = calloc(lineK_num, sizeof(FILE *));
    for (int i = 0; i < lineK_num; i++)
        res_file[i] = fopen(write_file_name, "wb");
    // src_array每一行对应的CELL个数(已对齐)
    int64_t read_cell_num = 0;
    // 保留未使用的冗余块需要的额外空间
    int64_t src_extra_cell_num = 0;
    // 根据冗余块推导原始块时，所需要的额外空间
    int64_t res_extra_cell_num = (lineK_num - 1) * (lineK_num / 2);
    /* 源文件的读取划分点，用于将src_array的每一行映射到文件对应位置，单位为byte */
    int64_t file_read_points[lineK_num];
    /* 记录每一个冗余块的偏移单位 */ 
    int64_t block_offests[lineK_num];
    /* 记录每一个冗余块的偏移总和 */
    int64_t block_extras[lineK_num];
    /* 冗余块的读取绝对偏移量 */ 
    int64_t src_read_offests[lineK_num];
    /* 原始块的写入相对偏移量 */ 
    int64_t res_write_offests[lineK_num];
    /* 原始块的写入绝对偏移量 */ 
    int64_t res_write_start[lineK_num];
    /* xorK异或起始点的绝对偏移量 */
    int64_t res_extras[lineK_num];
    /* 避免无效块写入结果数组 */ 
    int64_t not_write_res_num[lineK_num];
    /* 有效计算块的数据量 */ 
    int64_t calculate_left_num[lineK_num];
    /* 冗余块需要读取的数据量 */
    int64_t read_left_num[lineK_num];
    // 表示中间的冗余块的位置
    int64_t the_middle = (lineK_num - 1) / 2;


    // 计算每个冗余块的偏移单位
    for (int i = 0; i < lineK_num; i++)
        block_offests[i] = i - the_middle;
    // 计算每个冗余块的偏移总和
    for (int i = 0; i < lineK_num; i++)
        block_extras[i] = (lineK_num - 1) * labs(block_offests[i]);

    
    // 计算写入相对偏移量
    res_write_offests[the_middle] = 0;
    for (int i = the_middle + 1; i < lineK_num; i++)
        res_write_offests[i] = res_write_offests[i - 1] + block_offests[i - 1];
    for (int i = 1; i <= the_middle; i++)
        res_write_offests[the_middle - i] = res_write_offests[the_middle + i];
    // 计算读取绝对偏移量
    for (int i = the_middle; i < lineK_num; i++)
        src_read_offests[i] = (lineK_num - i - 1) * block_offests[i] + res_write_offests[i];
    for (int i = 0; i < the_middle; i++)
        src_read_offests[i] = i * labs(block_offests[i]) + res_write_offests[i];
    src_extra_cell_num = src_read_offests[lineK_num - 1];


    // 计算写入绝对偏移量
    for (int i = 0; i < lineK_num; i++)
        res_write_start[i] = res_extra_cell_num - src_extra_cell_num + res_write_offests[i];
    // 计算异或起始点的绝对偏移量
    for (int i = 0; i < lineK_num; i++)
        res_extras[i] = res_write_start[i] - i * block_offests[i];
    // 计算每行无效块个数
    for (int i = 0; i < lineK_num; i++)
        not_write_res_num[i] = src_extra_cell_num - res_write_offests[i];


    // 计算真实文件大小
    for (int i = 0; i < lineK_num; i++)
        true_cell_num -= block_extras[i];
    // 文件划分后每一行的CELL个数，包含了补零行
    int64_t line_cell_num = true_cell_num / lineK_num;
    // 计算有效计算块的数据量
    for (int i = 0; i < lineK_num; i++)
        calculate_left_num[i] = line_cell_num;
    // 计算需要读取的冗余块数据量
    for (int i = the_middle; i < lineK_num; i++)
        read_left_num[i] = (lineK_num - i - 1) * block_offests[i] + line_cell_num;
    for (int i = 0; i < the_middle; i++)
        read_left_num[i] = i * labs(block_offests[i]) + line_cell_num;
    // 源文件的读取划分点初始化，每个冗余块的大小不同
    file_read_points[0] = 0;
    for (int i = 1; i < lineK_num; i++)
        file_read_points[i] = file_read_points[i - 1] + (line_cell_num + block_extras[i - 1]) * BYTE_PERCELL;
    // 结果文件的写入划分点初始化，用于写入时确定每个写入点位置
    for (int64_t i = 0; i < lineK_num; i++)
        fseek(res_file[i], i * line_cell_num * BYTE_PERCELL, SEEK_SET);


    // 文件过大的情况，Memory装不下整个文件，需要多次读取文件内容
    // 将MAXMEMORY对齐到CELLSIZE，再将对应的CELL个数对齐到K，再计算src_array每行的CELL个数
    if (2 * all_cell_num * CELLSIZE > MAXMEMORY)
        read_cell_num = UPTO_K(UPTO_K(MAXMEMORY, CELLSIZE) / CELLSIZE, lineK_num) / lineK_num;
    // 文件不超限的情况，Memory装的下整个文件，无需多次读取文件内容
    // src_array每行的CELL个数就等于文件划分后每一行的CELL个数
    else
        read_cell_num = line_cell_num;
    // 需要写入的数据量
    int64_t write_left_size = line_cell_num;
    // 计算出原始块需要的次数
    int64_t calculate_times = UPTO_K(line_cell_num + src_extra_cell_num, read_cell_num) / read_cell_num;
    // 原始块起始点的绝对位置
    int64_t write_start = res_extra_cell_num;
    // 零偏有效数据终点的绝对位置
    int64_t write_over = res_extra_cell_num - src_extra_cell_num + read_cell_num;
    

    // src_array源数组，res_array结果数组的内存初始化
    // src_array的结构为(lineK_num行，read_cell_num列)的二维数组
    src_array = (__uint128_t **)aligned_alloc(lineK_num * sizeof(__uint128_t *), CELLSIZE);
    for (int i = 0; i < lineK_num; i++)
        src_array[i] = (__uint128_t *)aligned_alloc((read_cell_num + src_extra_cell_num) * sizeof(__uint128_t), CELLSIZE);
    // res_array的结构为(lineK_num行，read_cell_num + extra_cell_num列)的二维数组
    res_array = (__uint128_t **)aligned_alloc(lineK_num * sizeof(__uint128_t), CELLSIZE);
    for (int i = 0; i < lineK_num; i++)
        res_array[i] = (__uint128_t *)aligned_alloc((read_cell_num + res_extra_cell_num) * sizeof(__uint128_t), CELLSIZE);
    for (int i = 0; i < lineK_num; i++)
        memset(src_array[i], 0, (read_cell_num + src_extra_cell_num) * BYTE_PERCELL);
    for (int i = 0; i < lineK_num; i++)
        memset(res_array[i], 0, (read_cell_num + res_extra_cell_num) * BYTE_PERCELL);

    time_counter time_c;
    time_counter_init(&time_c);
    time_counter_begin(&time_c);
    for (int i = 0; i < calculate_times; i++)
    {
        // 此for循环进行数据读取操作
        for (int j = 0;  j < lineK_num; j++)
        {
            // 每次读取的数据量不能超过read_cell_num或剩余数据line_cell_num - read_counter
            int64_t need_read_size = read_cell_num < read_left_num[j] ? read_cell_num : read_left_num[j];
            fseek(src_file, file_read_points[j], SEEK_SET);
            fread(src_array[j] + src_extra_cell_num, sizeof(__uint128_t), need_read_size, src_file);
            file_read_points[j] += need_read_size * BYTE_PERCELL;
            read_left_num[j] -= need_read_size;
        }

        // 此for循环进行原始块的计算，j是第j列，k是第k块
        for (int j = 0; j < read_cell_num; j++)
        {
            for (int k = 0; k < the_middle; k++)
                if (not_write_res_num[k] > 0)
                    not_write_res_num[k]--;
                else if(calculate_left_num[k] > 0)
                {
                    calculate_left_num[k]--;
                    for (int64_t tt = 1; tt < TEST_TIMES; tt++)
                        simd_based_xorK(res_array, lineK_num, res_extras[k] + j, block_offests[k], src_array[k][src_read_offests[k] + j]);
                    res_array[k][res_write_start[k] + j] = simd_based_xorK(res_array, lineK_num, 
                        res_extras[k] + j, block_offests[k], src_array[k][src_read_offests[k] + j]);
                }
            for (int k = lineK_num - 1; k > the_middle ; k--)
                if (not_write_res_num[k] > 0)
                    not_write_res_num[k]--;
                else if(calculate_left_num[k] > 0)
                {
                    calculate_left_num[k]--;
                    for (int64_t tt = 1; tt < TEST_TIMES; tt++)
                        simd_based_xorK(res_array, lineK_num, res_extras[k] + j, block_offests[k], src_array[k][src_read_offests[k] + j]);
                    res_array[k][res_write_start[k] + j] = simd_based_xorK(res_array, lineK_num, 
                        res_extras[k] + j, block_offests[k], src_array[k][src_read_offests[k] + j]);
                }
            if (not_write_res_num[the_middle] > 0)
                not_write_res_num[the_middle]--;
            else if(calculate_left_num[the_middle] > 0)
            {
                calculate_left_num[the_middle]--;
                for (int64_t tt = 1; tt < TEST_TIMES; tt++)
                    simd_based_xorK(res_array, lineK_num, res_extras[the_middle] + j, block_offests[the_middle], src_array[the_middle][src_read_offests[the_middle] + j]);
                res_array[the_middle][res_write_start[the_middle] + j] = simd_based_xorK(res_array, lineK_num, 
                    res_extras[the_middle] + j, block_offests[the_middle], 
                    src_array[the_middle][src_read_offests[the_middle] + j]);
            }
        }
        
        // 进行数据块的写入
        int64_t need_write_size = 0 > write_over - write_start ? 0 : write_over - write_start;
        need_write_size = write_left_size < need_write_size ? write_left_size : need_write_size;
        for (int j = 0; j < lineK_num; j++)
            fwrite(res_array[j] + write_start, sizeof(__uint128_t), need_write_size, res_file[j]);
        write_left_size -= need_write_size;
        // 此处进行extra_cell的更新
        shiftarray(res_array, lineK_num, res_extra_cell_num, read_cell_num);
        shiftarray(src_array, lineK_num, src_extra_cell_num, read_cell_num);
        // 更新原始块起始点相对位置
        write_start = write_start > write_over ? 
            write_start - read_cell_num : res_extra_cell_num - src_extra_cell_num;
        for (int i = 0; i < lineK_num; i++)
            memset(src_array[i] + src_extra_cell_num, 0, read_cell_num * BYTE_PERCELL);
        for (int i = 0; i < lineK_num; i++)
            memset(res_array[i] + res_extra_cell_num, 0, read_cell_num * BYTE_PERCELL);
    }
    time_counter_end(&time_c);
    // printf("decode ver2, filename = %s,  k = %ld,  cost time = %lu\n", read_file_name, lineK_num, finish_time-start_time);
    dspf("./ente.txt", 4, file_size, lineK_num, get_run_time(time_c), true_cell_num, TEST_TIMES);


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