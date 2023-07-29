# Kernighan-Lin
Multi-module division based onKernighan-Lin algorithm

参考[陈仪香](https://faculty.ecnu.edu.cn/_s43/cyx/main.psp)老师课程：软硬件协同设计--第7章内容多模块划分方法 

## 任务图数据生成

参照  [TGFF](https://robertdick.org/projects/tgff/index.html) 进行任务图生成。

### a. 用于生成文件的tgffopt文件：（文件名 t22.tgffopt）

```
tg_cnt 5
task_cnt 20 5
task_degree 4 3
period_laxity 1
period_mul 1
tg_write
eps_write
vcg_write

table_label COMMUN
table_cnt 5
table_attrib price 80 20
type_attrib exec_time 50 20
trans_write
```

之后对应目录运行命令：tgff filename

<img src="figure\image-20230729200737570.png" alt="image-20230729200737570" style="zoom:80%;" />

filename 不需要添加后缀

之后会得到同名3个文件： .eps  .tgff  .vcg

**t22.eps** 如下，选取任务个数为21的任务图（配置生成相关，生成的是一定范围内的任务图）

<img src="figure\image-20230729200934530.png" alt="image-20230729200934530" style="zoom:80%;" />

数据需要使用 .tgff 文件内数据，为了是其符合本算法数据输入情况，对其进行预处理。

### b. 使用 python 脚本进行数据预处理（t22.tgff  ==> t2.csv)

python 脚本文件如下：(脚本名：**t.py**)

```
#!/usr/bin/python
#argv[1] is a filename; argv[2] is the number of graph
import re,sys,string
fp = open(sys.argv[1], "r")
all = fp.read()
fp.close
for i in range(int(sys.argv[2])):
        pattern2 = re.compile(r'FROM\s*t' + str(i) + '_(\S*)\s*TO\s*t' + str(i) + '_(\S*)\s*TYPE\s*(\S*)')
        fp_o = open('t' + str(i) + '.csv', "w")
        # fp_o.write('Source,Target,price\n')
        for match2 in pattern2.finditer(all):
                print (match2.group(1,2,3))
                fp_o.write(match2.group(1) + ' ' + match2.group(2) + ' ' + match2.group(3) + '\n')
        fp_o.close
```

对应文件目录运行： 

<img src="figure\image-20230729202043180.png" alt="image-20230729202043180" style="zoom:80%;" />

t22.tgff 和 5 为脚本的传入参数，分别为要预处理的数据及处理图的数量。下面三元组分别为头任务、尾任务和任务间通信代价

得到符合算法需求的数据文件 **t1.csv**

```
0 1 18
0 2 12
1 2 13
2 3 40
0 3 8
1 3 27
3 4 9
4 5 8
1 5 40
3 5 44
2 5 0
2 6 30
4 7 15
3 7 2
6 7 23
7 8 41
7 9 21
8 10 39
7 11 2
9 11 35
10 11 48
5 11 28
6 12 14
4 12 28
8 12 40
11 13 30
9 13 42
13 14 32
13 15 18
12 16 10
12 17 18
16 18 25
16 19 33
16 20 27
```

同理可获得任务数分别为 21、 33、51的数据文件，对应任务图如下：

<img src="figure\21" alt="image-20230729204033136" style="zoom:67%;" />

<img src="figure\33" alt="image-20230729204736037" style="zoom:67%;" />



<img src="figure\51" alt="image-20230729205642478" style="zoom:67%;" />
