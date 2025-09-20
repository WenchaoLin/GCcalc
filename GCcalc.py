#!/usr/bin/env python3
"""
GC含量和GC偏斜计算工具

该程序用于计算FASTA格式序列文件中每个滑动窗口的GC含量和GC偏斜值。
支持自定义窗口大小和步长参数。

作者: Wenchao Lin
版权: Copyright 2025, Uniteomics Biotech Corp.
许可证: GPL
版本: 1.0.1
维护者: Wenchao Lin
邮箱: linwenchao@uniteomics.com
"""

from __future__ import division
import argparse
import sys
import os
from typing import Tuple, Iterator, TextIO

try:
    from Bio import SeqIO
except ImportError:
    print("错误: 需要安装Biopython库")
    print("请运行: pip install biopython")
    sys.exit(1)


# 常量定义
DEFAULT_WINDOW_SIZE = 1000
DEFAULT_STEP_SIZE = 1000
GC_BASES = {'G', 'C', 'g', 'c', 'S', 's'}  # S代表G或C的简并碱基
G_BASES = {'G', 'g'}
C_BASES = {'C', 'c'}
PRECISION = 4


class GCCalculator:
    """GC含量和GC偏斜计算器"""
    
    def __init__(self, window_size: int = DEFAULT_WINDOW_SIZE, 
                 step_size: int = DEFAULT_STEP_SIZE):
        """
        初始化计算器
        
        Args:
            window_size: 窗口大小
            step_size: 步长大小
        """
        self.window_size = window_size
        self.step_size = step_size
        self._validate_parameters()
    
    def _validate_parameters(self) -> None:
        """验证参数的有效性"""
        if self.window_size <= 0:
            raise ValueError("窗口大小必须大于0")
        if self.step_size <= 0:
            raise ValueError("步长必须大于0")
        if self.window_size < self.step_size:
            print("警告: 窗口大小小于步长，可能存在重叠")
    
    def calculate_gc_content(self, sequence: str) -> float:
        """
        计算序列的GC含量
        
        Args:
            sequence: 输入序列
            
        Returns:
            GC含量 (0-1之间的浮点数)
        """
        if not sequence:
            return 0.0
        
        sequence_upper = sequence.upper()
        gc_count = sum(1 for base in sequence_upper if base in GC_BASES)
        return round(gc_count / len(sequence), PRECISION)
    
    def calculate_gc_skew(self, sequence: str) -> float:
        """
        计算序列的GC偏斜
        
        Args:
            sequence: 输入序列
            
        Returns:
            GC偏斜值 (-1到1之间的浮点数)
        """
        if not sequence:
            return 0.0
        
        sequence_upper = sequence.upper()
        g_count = sum(1 for base in sequence_upper if base in G_BASES)
        c_count = sum(1 for base in sequence_upper if base in C_BASES)
        
        total_gc = g_count + c_count
        if total_gc == 0:
            return 0.0
        
        skew = (g_count - c_count) / total_gc
        return round(skew, PRECISION)
    
    def process_sequence(self, sequence: str, sequence_id: str) -> Iterator[Tuple[str, int, int, float, float]]:
        """
        处理单个序列，生成滑动窗口结果
        
        Args:
            sequence: 序列字符串
            sequence_id: 序列标识符
            
        Yields:
            (序列ID, 起始位置, 结束位置, GC含量, GC偏斜) 的元组
        """
        seq_length = len(sequence)
        
        for start_pos in range(0, seq_length, self.step_size):
            end_pos = min(start_pos + self.window_size, seq_length)
            window_seq = sequence[start_pos:end_pos]
            
            # 确保窗口大小足够
            if len(window_seq) < self.window_size:
                # 对于最后一个不完整的窗口，可以选择跳过或处理
                # 这里选择处理，保持与原程序行为一致
                pass
            
            gc_content = self.calculate_gc_content(window_seq)
            gc_skew = self.calculate_gc_skew(window_seq)
            
            # 位置从1开始计数（生物信息学惯例）
            yield (sequence_id, start_pos + 1, end_pos, gc_content, gc_skew)


def validate_fasta_file(filepath: str) -> None:
    """
    验证FASTA文件的有效性
    
    Args:
        filepath: 文件路径
        
    Raises:
        FileNotFoundError: 文件不存在
        ValueError: 文件格式无效
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"文件不存在: {filepath}")
    
    if not os.path.isfile(filepath):
        raise ValueError(f"路径不是文件: {filepath}")
    
    # 尝试读取第一个序列来验证格式
    try:
        with open(filepath, 'r') as handle:
            first_record = next(SeqIO.parse(handle, 'fasta'), None)
            if first_record is None:
                raise ValueError("FASTA文件中没有找到序列")
    except Exception as e:
        raise ValueError(f"无效的FASTA文件格式: {e}")


def process_fasta_file(filepath: str, calculator: GCCalculator, 
                      output_file: TextIO = None) -> None:
    """
    处理FASTA文件并输出结果
    
    Args:
        filepath: FASTA文件路径
        calculator: GC计算器实例
        output_file: 输出文件句柄，None表示输出到标准输出
    """
    try:
        with open(filepath, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequence_id = record.id
                sequence = str(record.seq)
                
                for result in calculator.process_sequence(sequence, sequence_id):
                    seq_id, start, end, gc_content, gc_skew = result
                    output_line = f"{seq_id}\t{start}\t{end}\t{gc_content}\t{gc_skew}"
                    
                    if output_file:
                        output_file.write(output_line + '\n')
                    else:
                        print(output_line)
    
    except Exception as e:
        print(f"处理文件时出错: {e}", file=sys.stderr)
        raise


def create_argument_parser() -> argparse.ArgumentParser:
    """创建命令行参数解析器"""
    parser = argparse.ArgumentParser(
        description="计算FASTA序列的GC含量和GC偏斜",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  %(prog)s -f input.fa                    # 使用默认参数
  %(prog)s -f input.fa -w 500 -s 250     # 自定义窗口和步长
  %(prog)s -f input.fa -o output.txt     # 输出到文件
        """
    )
    
    parser.add_argument(
        '-f', '--file',
        dest='filename',
        required=True,
        help='输入的FASTA格式文件'
    )
    
    parser.add_argument(
        '-w', '--window',
        dest='window_size',
        type=int,
        default=DEFAULT_WINDOW_SIZE,
        help=f'窗口大小 (默认: {DEFAULT_WINDOW_SIZE})'
    )
    
    parser.add_argument(
        '-s', '--step',
        dest='step_size',
        type=int,
        default=DEFAULT_STEP_SIZE,
        help=f'步长大小 (默认: {DEFAULT_STEP_SIZE})'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='output_file',
        help='输出文件路径 (默认输出到标准输出)'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 2.0.0'
    )
    
    return parser


def main() -> None:
    """主函数"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    try:
        # 验证输入文件
        validate_fasta_file(args.filename)
        
        # 创建计算器
        calculator = GCCalculator(args.window_size, args.step_size)
        
        # 处理文件
        if args.output_file:
            with open(args.output_file, 'w') as output_handle:
                process_fasta_file(args.filename, calculator, output_handle)
            print(f"结果已保存到: {args.output_file}", file=sys.stderr)
        else:
            process_fasta_file(args.filename, calculator)
    
    except (FileNotFoundError, ValueError) as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n程序被用户中断", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"未预期的错误: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
