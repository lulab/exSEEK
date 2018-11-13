# zhuoer's part of the exSeek project


# useful command

```bash
lulab
cd /BioII/lulab_b/dongzhuoer/
cd /home/chenxupeng/projects/exseek/

scp -r -P 10358    zhuoer@166.111.156.58:/home/chenxupeng/projects/exseek/ .
```

# useful location


# to xvpeng

```r
m <- read.csv(matrix_path, sep='\t',row.names = 1, header=TRUE)
```

- `read.csv(sep='\t')` 不如直接用 `read.table()`。
- 只读前 1000 行，方便后面计算

```r
m <- read.table(matrix_path,row.names = 1, header=TRUE, nrows = 10000)
```