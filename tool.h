#include <emmintrin.h>
#include <stdalign.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <sys/time.h>
// #include <windows.h>

#define BYTESIZE (8LL)
#define CELLSIZE (128LL)
#define BYTE_PERCELL (CELLSIZE / BYTESIZE)
#define MAXMEMORY (4LL * 1024 * 1024 * 1024 * BYTESIZE) // 4GB，物理存储内存
#define MAXCACHE (256LL * 1024 * BYTESIZE) // 256KB，物理高速缓存
#define UPTO_K(num, K) (((num + (K) - 1LL) / (K)) * (K))
#define MICRO (1000000LL) // 秒转微秒
#define TEST_SIZE (32LL * 1024 * 1024 * BYTESIZE)
#define TEST_TIMES (TEST_SIZE / file_size)

typedef struct cal_args
{
    int64_t file_size;
    int64_t line_read_times;
    int64_t read_cell_num;
    int64_t lineK_num;
    int64_t index;
    int64_t *file_write_nums;
    int64_t *res_extras;
    int64_t *block_offests;
    alignas(CELLSIZE) __uint128_t **src_array;
    alignas(CELLSIZE) __uint128_t **res_array;
    sem_t *read_done;
    sem_t *cal_done;
    sem_t *write_done;
}cal_args;

typedef struct rw_args
{
    int64_t line_read_times;
    int64_t read_cell_num;
    int64_t line_cell_num;
    int64_t extra_cell_num;
    int64_t lineK_num;
    int64_t *file_read_points;
    int64_t *file_write_nums;
    alignas(CELLSIZE) __uint128_t **src_array;
    alignas(CELLSIZE) __uint128_t **res_array;
    FILE *src_file; 
    FILE **res_file;
    sem_t *read_done;
    sem_t *cal_done;
    sem_t *write_done;
}rw_args;

typedef struct time_counter
{
    struct timeval begin_time;
    struct timeval end_time;
    // LARGE_INTEGER freq;
}time_counter;

// 初始化rw参数
void rw_args_init(rw_args *irw_args, int64_t line_read_times, int64_t read_cell_num, int64_t line_cell_num, int64_t extra_cell_num, 
    int64_t lineK_num, int64_t *file_read_points, int64_t *file_write_nums, __uint128_t **src_array, __uint128_t **res_array, 
    FILE *src_file, FILE **res_file, sem_t *read_done, sem_t *cal_done, sem_t *write_done);

// 初始化cal参数
void cal_args_init(cal_args *ical_args, int64_t file_size, int64_t line_read_times, int64_t read_cell_num, int64_t lineK_num, int64_t index, 
    int64_t *file_write_nums, int64_t *res_extras, int64_t *block_offests, __uint128_t **src_array, __uint128_t **res_array, 
    sem_t *read_done, sem_t *cal_done, sem_t *write_done);

// 初始化计时器
void time_counter_init(time_counter *time_c);

// 计时器开始计时
void time_counter_begin(time_counter *time_c);

// 计时器结束计时
void time_counter_end(time_counter *time_c);

// 获取运行时间，微秒级精度
double get_run_time(time_counter time_c);

// 返回系统时间
// uint64_t rdtscp();

// 对齐的malloc，注意未初始化
void *aligned_alloc(size_t _Size, size_t _Alignment);

// 对齐的free
void aligned_free(void *_Memory);

// 打印从target指向的位置开始，length长(单位为BYTE)的数据的16进制表示
// 参数1：target待打印的数据的起始地址
// 参数2：length待打印的数据的长度
void printbit(void *target, int length);

// 将arg二维数组(lineK_num行)的最后extra_cell_num列移到起始列
// 参数1：arg待移位的数组的起始地址
// 参数2：lineK_num数组的行数
// 参数3：shift_cell_num数组要移动的列数
// 参数4：read_cell_num数组的预读取数据列数
void shiftarray(__uint128_t **arg, long lineK_num, long shift_cell_num, long read_cell_num);

// line_num个128 bit的二进制数进行异或(xor)操作
// 选取arg数组的(i行, line_shift + rowshift * i列)的元素
// 将选取的元素与res进行异或，i从0到line_num - 1
// 参数1：arg待异或的数组的起始地址
// 参数2：lineK_num数组的行数
// 参数3：line_shift固定行偏移量
// 参数4：row_shift列增行偏移量
__uint128_t xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift);

// xorK的偏置版本
// 参数1：arg待异或的数组的起始地址
// 参数2：lineK_num数组的行数
// 参数3：line_shift固定行偏移量
// 参数4：row_shift列增行偏移量
// 参数5：based初始偏置值
__uint128_t based_xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift, __uint128_t based);

// xorK的SIMD版本(每次计算128 bit)
// 参数1：arg待异或的数组的起始地址
// 参数2：lineK_num数组的行数
// 参数3：line_shift固定行偏移量
// 参数4：row_shift列增行偏移量
__uint128_t simd_xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift);

// xorK的偏置SIMD版本(每次计算128 bit)
// 参数1：arg待异或的数组的起始地址
// 参数2：lineK_num数组的行数
// 参数3：line_shift固定行偏移量
// 参数4：row_shift列增行偏移量
// 参数5：based初始偏置值
__uint128_t simd_based_xorK(__uint128_t **arg, long lineK_num, int line_shift, int row_shift, __uint128_t based);

// 多线程执行异或计算部分
// 参数1：arg参数数组
void *cal_parallel(void *arg);

// 执行数据读取写入部分
// 参数1：arg参数数组
void *rw_parallel(void *arg);

// 去除.backup.xx文件扩展名
char *realname(char *name);

// 获取encode时使用的k
long getK(char *name);

// 使用k=3的矩阵进行直接运算，不考虑文件大小等问题
// 参数1：line_cell_num每行的元素个数
// 参数2：src_file源文件
// 参数2：res_file目标文件
void encodek3(long line_cell_num, FILE *src_file, FILE *res_file);

// 使用k=4的矩阵进行直接运算，不考虑文件大小等问题
// 参数1：line_cell_num每行的元素个数
// 参数2：src_file源文件
// 参数2：res_file目标文件
void encodek4(long line_cell_num, FILE *src_file, FILE *res_file);

// 使用k=5的矩阵进行直接运算，不考虑文件大小等问题
// 参数1：line_cell_num每行的元素个数
// 参数2：src_file源文件
// 参数2：res_file目标文件
void encodek5(long line_cell_num, FILE *src_file, FILE *res_file);

// 使用k=6的矩阵进行直接运算，不考虑文件大小等问题
// 参数1：line_cell_num每行的元素个数
// 参数2：src_file源文件
// 参数2：res_file目标文件
void encodek6(long line_cell_num, FILE *src_file, FILE *res_file);

// 编码消耗
// fsz bit单位，sum CELL单位
void espf(char *fna, int64_t ver, int64_t fsz, int64_t k, double ct, int64_t sum, int64_t test_times);

// 解码消耗
// fsz bit单位，sum CELL单位
void dspf(char *fna, int64_t ver, int64_t fsz, int64_t k, double ct, int64_t sum, int64_t test_times);