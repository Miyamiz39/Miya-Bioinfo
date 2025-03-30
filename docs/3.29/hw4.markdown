# 作业4

## 1. Bowtie 与 BWT（Burrows-Wheeler Transform）

*BWT变换的基本过程是：*
  
 *1. “旋转” (rotate): 每次将字符串向右旋转一个字符（将末尾的字符移动到字符串的开头），得到一个字符串数组。*

 *2. 排序：对字符串数组的第一列按照字典序排序。若相同，再按照下一列，依次类推。*

 *3. 拼接：将排序后的字符串数组的<u>最后一列</u>拼接起来，得到结果。*

### BWT 提高运算速度的性质

- **压缩数据**：BWT变换得到的结果，使得相同字符聚集。例如，` banana `
经过变换后成为` bnnaaa `，游程编码可直接压缩为` "b1n2a3" `

- **无损压缩**：BWT变换的逆变换是BWT变换，即BWT(BWT(x)) = x。

    *使BWT(T)可逆的性质是 **LF映射** (Last-First Mapping)：
在最后一列第i次出现的字符，就是第一列中的第i次出现的那个字符。*

- Bowtie利用这一性质，将查询序列（read）从右向左逐字符比对（称为反向搜索）。


### Bowtie 节省内存的策略

1. 在匹配到序列后，为了快速定位序列位置，使用了 **里程碑（Milestones）**索引，即每一定的距离，添加一个索引。同时结合**回溯搜索**，从当前字符回溯至前一个索引，并计算步数。

2. 在为BWT转换后的序列建立索引（为了确定字符的顺位）时，使用了类似的**检查点（Checkpoints）**策略。每隔一定距离，记录这个字符的顺位。再结合**回溯搜索**。

    **里程碑**和**检查点**使得Bowtie算法的总内存消耗仅为 1.65 倍于基因组序列。其中：
    |||
    |:---:|:---:|
    | BWT | 1 x |
    |Milestones|~0.5 x|
    |Checkpoints|~0.15 x|
    |||
---

## 2. 

```bash 
bowtie -v 2 -m 10 --best --strata BowtieIndex/YeastGenome -f THA2.fa -S THA2.sam

#返回结果：
# reads processed: 1250
# reads with at least one reported alignment: 1158 (92.64%)
# reads that failed to align: 77 (6.16%)
# reads with alignments suppressed due to -m: 15 (1.20%)
```

- [比对结果 THA2.sam](https://raw.githubusercontent.com/Miyamiz39/Miya-Bioinfo/main/docs/3.29/THA2.sam)

- 各染色体上的 mapping 数量

    ```
     92 *           #未mapping的reads
     18 chrI
     51 chrII
     15 chrIII
  194 chrIV
     25 chrIX
     12 chrmt
     33 chrV
     17 chrVI
    125 chrVII
     68 chrVIII
     71 chrX
     56 chrXI
    169 chrXII
     67 chrXIII
     58 chrXIV
    101 chrXV
     78 chrXVI
     ```
## 3. <sup>[1]</sup>
### 3.1 什么是sam/bam文件中的"CIGAR string"? 

 `CIGAR`位于比对部分 (alignment segment) 的第六列，数据类型是字符串。`CIGAR` 表示了reads的每一个碱基与genome的比对情况：

| 字符 | 含义 |
|:---:|:---:|
| M | 对齐的匹配|
| I | 插入
| D | 删除
| N | 跳过（跳过一个genome片段）
| S | soft-clipped
| H | hard-clipped
| P | padding
| = | 匹配
| X | 错配

*示例：*
```
Read:    ATG--CGT  (比对到参考序列)
Ref:     ATGTTTGT
CIGAR:   3M2D3M   (3匹配 + 2删除 + 3匹配)
```

### 3.2 "soft clip"的含义是什么？

在 CIGAR 字符串中以 `S` 表示。是由于**测序质量低、被遭到污染**等原因而被比对程序认为无法比对的碱基。这些碱基不影响 POS 位置，但仍会被保留在 `SEQ` 列之中。

相比之下，hard clip 是 CIGAR 字符串中以 `H` 表示的，它一般出现在比对结果中，表示测序接头等无法用于比对的碱基，不会被保留在 `SEQ` 中。

*示例：*


![](https://raw.githubusercontent.com/Miyamiz39/Miya-Bioinfo/main/docs/3.29/ali.png)

结果为

![](https://raw.githubusercontent.com/Miyamiz39/Miya-Bioinfo/main/docs/3.29/result.png)

### 3.3 什么是reads的mapping quality? 

MAPQ 位于文件的第 5 列，用 Phred 分数表示，定义如下：

```
MAPQ = -10 * log10(P) 
# P 为比对出错的概率
```

一般取值在 0 - 60 之间，越大表示比对质量越高。唯一比对的读段通常有**高** MAPQ（如 ≥30）。重复序列区域或低复杂度区域的读段往往有**低** MAPQ（如 ≤10）。

### 3.4 通过比对结果推断参考序列？

若在比对结果最后一列中记录了 `MD` tag，则能够推断出比对的参考序列。

`MD` 的数据类型是字符串。格式为：
```
[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*

# [0-9]+ :match的碱基数量
# [A-Z]  :单个错配的碱基
# \^[A-Z]+ : 参考中被删除的碱基
```

也就是说，`MD` 字段记录了 reads 和被匹配到的参考序列的差异情况。因此，结合 `CIGAR` 和 `SEQ` 就可以推断出参考序列。



## 4.
- [比对结果 THA2-bwa.sam](https://raw.githubusercontent.com/Miyamiz39/Miya-Bioinfo/main/docs/3.29/THA2-bwa.sam)


## 5. Genome Browser的使用

可视化比对结果 `THA2.sam`

![](https://raw.githubusercontent.com/Miyamiz39/Miya-Bioinfo/main/docs/3.29/igv_snapshot.png)


---

### 参考资料：

[1] https://github.com/samtools/hts-specs