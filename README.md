## 编码
gcc encode3.c -o en3 -pthread
./en3 {filename} {k}

## 解码
gcc decode4.c -o de4 -pthread
./de4 {filename} {k}
