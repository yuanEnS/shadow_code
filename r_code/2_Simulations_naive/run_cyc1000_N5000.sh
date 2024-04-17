#!/bin/bash

### 将本次作业计费到导师课题组，tutor_project改为导师创建的课题组名
#SBATCH --comment=sos_weili

### 给你这个作业起个名字，方便识别不同的作业
#SBATCH --job-name=200-times-repeat-N-8000
  
### 指定该作业需要多少个节点
### 注意！没有使用多机并行（MPI/NCCL等），下面参数写1！不要多写，多写了也不会加速程序！
#SBATCH --nodes=1

### 指定该作业需要多少个CPU核心
### 注意！一般根据队列的CPU核心数填写，比如cpu队列64核，这里申请<=64核！
#SBATCH --ntasks=64

### 指定该作业在哪个队列上执行
### 目前可用的CPU队列有 cpu/fat
### cpu队列有64核，fat队列有128核
#SBATCH --partition=cpu64c

### 以上参数用来申请所需资源
### 以下命令将在计算节点执行

### 加载Anaconda
export PATH=/opt/app/anaconda3/bin:$PATH

### 激活R环境
source activate shadow

### 执行你的作业
Rscript Repeat_cyc1000_N5000_sace.R
