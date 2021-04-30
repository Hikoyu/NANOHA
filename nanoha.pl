#!/usr/bin/env perl

# Copyright (c) 2020-2021 Hikoyu Suzuki
# This software is released under the MIT License.

use strict;
use warnings;
use Getopt::Std;
use threads;

# ソフトウェアを定義
### 編集範囲 開始 ###
my $software = "nanoha.pl";	# ソフトウェアの名前
my $version = "ver.2.0.0";	# ソフトウェアのバージョン
my $note = "NANOHA is Network-based Assortment of Noisy On-target reads for High-accuracy Alignments.\n  This software assorts on-target PacBio/Nanopore reads such as target amplicon sequences.";	# ソフトウェアの説明
my $usage = "<required items> [optional items]";	# ソフトウェアの使用法 (コマンド非使用ソフトウェアの時に有効)
### 編集範囲 終了 ###

# コマンドを定義
my %command;
### 編集範囲 開始 ###
$command{"sketch"} = "Sketch out sequence reads";
$command{"build"} = "\tBuild sequence similarity graph";
$command{"assort"} = "Assort sequence reads based on sequence similarity graph";
$command{"unify"} = "\tUnify sequence reads in the same cluster";
$command{"convert"} = "Convert sequence reads from FASTA format to FASTQ format";
$command{"dump"} = "\tDump k-mer minimizers counts";

# コマンドを追加
### 編集範囲 終了 ###
my @command_list = sort(keys(%command));

# 指定されたコマンドを確認
my $specified_command = shift(@ARGV) if @command_list and @ARGV;
&exception::error("unknown command: $specified_command") if defined($specified_command) and !exists($command{$specified_command});

# 共通オプションを定義
my %option;
### 編集範囲 開始 ###
# オプションを追加
### 編集範囲 終了 ###

# コマンドごとのオプション定義を取得
&{\&{"${specified_command}::define"}} if $specified_command;
my @option_list = sort(keys(%option));

# ヘルプを表示 (引数未指定時)
&exception::help unless @ARGV or -p STDIN;

# オプションの入力処理
my %opt;
$_ = join("", @option_list);
$_ =~ s/\s+\S+\s+/:/g;
getopts($_, \%opt);

# 未指定オプションのデフォルト値を入力
foreach (@option_list) {
	$opt{substr($_, 0, 1)} = substr($option{$_}, index($option{$_}, "[") + 1, index($option{$_}, "]") - index($option{$_}, "[") - 1) if $option{$_} =~ /\[.+\]$/ and !defined($opt{substr($_, 0, 1)});
}

### 編集範囲 開始 ###
# 追加のモジュールを宣言
use threads::shared;
use Thread::Queue 3.07;
use List::Util 1.54;
use List::MoreUtils 0.420_001;
use Math::BigFloat;
use Math::Random::MT 1.17;
use Number::AnyBase 1.60000;
use Scalar::Util;
use if $ENV{"NANOHA_XS"}, Inline => (C => Config => CC => exists($ENV{"CC"}) ? $ENV{"CC"} : 'cc', CCFLAGS => '-std=c99', DIRECTORY => $ENV{"NANOHA_INLINE_DIR"});
use if $ENV{"NANOHA_XS"}, Inline => (C => 'DATA', NAME => 'NANOHA::BUILD::LLCS::SISD');
use if $ENV{"NANOHA_XS"}, Inline => (C => 'DATA', NAME => 'NANOHA::BUILD::LLCS::SIMD', CCFLAGS => ['-msse41', '-mavx2', '-mavx512f -mavx512bw']->[scalar(@{[`sysctl machdep.cpu.leaf7_features` =~ /^(?=.*AVX512F)(?=.*AVX512BW)|AVX2/g]})], AUTO_INCLUDE => ['#include <x86intrin.h>', '#define word_size_index ' . [1, 2, 3]->[scalar(@{[`sysctl machdep.cpu.leaf7_features` =~ /^(?=.*AVX512F)(?=.*AVX512BW)|AVX2/g]})]]);
use if $ENV{"NANOHA_XS"}, Inline => (CPP => Config => CC => exists($ENV{"CXX"}) ? $ENV{"CXX"} : 'c++', CCFLAGS => '-std=c++14 -march=native', DIRECTORY => $ENV{"NANOHA_INLINE_DIR"});
use if $ENV{"NANOHA_XS"}, Inline => (CPP => 'DATA', NAME => 'NANOHA::SKETCH', AUTO_INCLUDE => ['#undef seed', '#include <vector>', '#include <unordered_map>', '#include <random>', '#include <algorithm>']);
use if $ENV{"NANOHA_XS"}, Inline => (CPP => 'DATA', NAME => 'NANOHA::BUILD::SORT', AUTO_INCLUDE => ['#undef seed', '#include <string>', '#include <unordered_map>']);
use if $ENV{"NANOHA_XS"}, Inline => (CPP => 'DATA', NAME => 'NANOHA::ASSORT', AUTO_INCLUDE => ['#undef seed', '#include <vector>', '#include <unordered_map>', '#include <random>', '#include <algorithm>', '#include <numeric>']);
use if $ENV{"NANOHA_XS"} && $ENV{"NANOHA_SPOA"}, Inline => (CPP => 'DATA', NAME => 'NANOHA::UNIFY', AUTO_INCLUDE => ['#include <x86intrin.h>', '#include <spoa/spoa.hpp>'], INC => exists($ENV{"SPOA_INC"}) ? "-I$ENV{SPOA_INC}" : '-I/usr/local/inculde', LIBS => [exists($ENV{"SPOA_LIB"}) ? "-Wl,-rpath,$ENV{SPOA_LIB} -L$ENV{SPOA_LIB}" : '-Wl,-rpath,/usr/local/lib -L/usr/local/lib', '-lspoa']);
no warnings 'portable';

# 定数を定義
use constant {
	max_num_reads => 4294967295,		# max_num_reads => リード数上限値
	max_depth => 65535,					# max_depth => 深度上限値
};

# 塩基配列の基数変換則を定義
my $convert = Number::AnyBase->new_dna;

# 処理を追加
### 編集範囲 終了 ###

# メインルーチンを実行
&main;
exit(0);

# メインルーチン
sub main {
	### 編集範囲 開始 ###
	# 指定された共通オプションを確認
	
	# 処理を追加
	### 編集範囲 終了 ###
	
	# コマンドの実行 (コマンド指定時)
	&{\&{"${specified_command}::body"}} if $specified_command;
	
	### 編集範囲 開始 ###
	# 処理を追加
	### 編集範囲 終了 ###
	return(1);
}

### 編集範囲 開始 ###
# サブルーチンを追加
### 編集範囲 終了 ###

## ここから例外処理のパッケージ ##
package exception;

# ヘルプ表示
sub help {
	print STDERR "$software ";
	print STDERR $specified_command ? $specified_command : $version;
	print STDERR "\n\nFunctions:\n  $note\n\nUsage:\n  $software ";
	if (!$specified_command and @command_list) {
		print STDERR "<command>\n";
		print STDERR "\nCommand:\n";
		foreach (@command_list) {print STDERR "  $_\t$command{$_}\n";}
	}
	else {
		print STDERR "$specified_command " if $specified_command;
		print STDERR "[options] " if @option_list;
		print STDERR "$usage\n";
		print STDERR "\nOptions:\n" if @option_list;
		foreach (@option_list) {print STDERR "  -$_\t$option{$_}\n";}
	}
	exit(0);
}

# エラー表示
sub error {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Error: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	threads->tid or map {$_->detach} threads->list;	# threadsモジュールを使用する場合はアンコメント
	exit(1);
}

# 注意表示
sub caution {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Caution: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	return(1);
}

### 編集範囲 開始 ###
## ここからsketchコマンドのパッケージ ##
package sketch;

# コマンドとオプションを定義
sub define {
	$note = "Sketch out sequence reads under specified conditions.";
	$usage = "<STDIN|in1.fq> [in2.fq ...]";
	$option{"o STR "} = "Output file prefix [nanoha]";
	$option{"p INT "} = "Number of parallel worker threads <1-> [1]";
	$option{"k INT "} = "Size of k-mer <5-15> [15]";
	$option{"n INT "} = "Number of k-mer minimizers to be generated from each sequence read <1-255> [10]";
	$option{"u INT "} = "Maximum amount of sequence reads to be loaded <1-" . main::max_num_reads . "> [" . main::max_num_reads . "]";
	$option{"5 INT "} = "Number of bases trimmed from 5' (left) end of each sequence read <0-> [0]";
	$option{"3 INT "} = "Number of bases trimmed from 3' (right) end of each sequence read <0-> [0]";
	$option{"s"} = "\tUse strand-specific sequence reads";
	$option{"w"} = "\tUse 2-byte line feed code (CR+LF) for input files";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT 7-15: -k $opt{k}") if $opt{"k"} !~ /^\d+$/ or $opt{"k"} < 7 or $opt{"k"} > 15;
	&exception::error("specify INT 1-255: -n $opt{n}") if $opt{"n"} !~ /^\d+$/ or $opt{"n"} < 1 or $opt{"n"} > 255;
	&exception::error("specify INT 1-" . main::max_num_reads . ": -u $opt{u}") if $opt{"u"} !~ /^\d+$/ or $opt{"u"} == 0 or $opt{"u"} > main::max_num_reads;
	&exception::error("specify INT >= 0: -5 $opt{5}") if $opt{"5"} !~ /^\d+$/;
	&exception::error("specify INT >= 0: -3 $opt{3}") if $opt{"3"} !~ /^\d+$/;
	&exception::caution("even k-mer size not recommended: -k $opt{k}") unless $opt{"k"} & 0x01;
	
	# 使用コードを表示
	print STDERR "Use ", $ENV{"NANOHA_XS"} ? "XS" : "Perl", " code\n";
	
	# 使用パラメータを表示
	print STDERR "Size of k-mer: $opt{k}\n";
	print STDERR "Number of k-mer minimizers: $opt{n}\n";
	print STDERR "Number of bases trimmed from 5'-end: $opt{5}\n";
	print STDERR "Number of bases trimmed from 3'-end: $opt{3}\n";
	
	# 入力ファイルを確認
	&exception::error("input file not specified") unless @ARGV or -p STDIN;
	&common::check_files(\@ARGV);
	
	# プロセスIDとプログラム開始時刻をファイルヘッダーに登録
	my $file_header = pack("NN", $$, $^T);
	
	# 最大minimizerカウント数を算出
	my $max_minimizer_count = 2 ** (32 - $opt{"k"} * 2);
	
	# 変数を宣言
	my @worker_threads = ();
	my $num_error_threads = 0;
	
	# 入出力キューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 入出力キューの要素数上限を定義
	$input->limit = $opt{"p"};
	$output->limit = $opt{"p"};
	
	# 指定されたワーカースレッド数で並列処理
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_threads[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューを読み込みながら処理
			while (defined(my $serialized_data = $input->dequeue)) {
				# 変数を宣言
				my %kmers = ();
				
				# リードIDとリードシーケンスを取得
				my ($read_id, $read_seq) = unpack("LA*", $serialized_data);
				
				# 順鎖のシーケンスに含まれるk-merを取得
				$kmers{$_}++ foreach map {substr($read_seq, $_, $opt{"k"})} 0..length($read_seq) - $opt{"k"};
				
				# 逆鎖のシーケンスに含まれるk-merを取得 (-s未指定時)
				$kmers{$_}++ foreach map {&common::complementary($_)} $opt{"s"} ? () : keys(%kmers);
				
				# weighted minimizerを生成
				my $weighted_minimizers = &generate_weighted_minimizers([List::Util::pairmap {$convert->to_dec("C" . $a) => $b} %kmers], $opt{"n"});
				
				# リードID及び各weighted minimizerをバイナリ形式で出力キューに追加
				$output->enqueue(pack("N*", $read_id, @{$weighted_minimizers}));
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# データストリームスレッドを作成
	## ここからデータストリームスレッドの処理 ##
	$worker_threads[$opt{"p"}] = threads::async {
		# このスレッドのみを終了可能に変更
		threads->set_thread_exit_only(1);
		
		# Nanoha Sequence Sketch (NSS) ファイルを作成
		open(NSS, ">", "$opt{o}.nss") or &exception::error("failed to make file: $opt{o}.nss");
		
		# NSSファイルをバイナリモードにする
		binmode(NSS);
		
		# ファイルヘッダーとk-merサイズ、ストランド特異性、メジャーバージョン及びminimizerの生成個数をバイナリ形式でNSSファイルに出力
		print NSS $file_header, pack("CC", [split(/\./, $version)]->[1] - 1 << 5 | !!$opt{"s"} << 4 | $opt{"k"}, $opt{"n"});
		
		# 変数を宣言
		my $seq_sketch_index = "";
		
		# 出力キューを読み込みながら処理
		while (defined(my $serialized_data = $output->dequeue)) {
			# NSSファイルのファイルポインタの位置をシーケンススケッチインデックスに登録
			vec($seq_sketch_index, vec($serialized_data, 0, 32), 64) = tell(NSS);
			
			# 各minimizerをバイナリ形式でNSSファイルに出力
			print NSS substr($serialized_data, 4);
		}
		
		# シーケンススケッチインデックス (先頭の8バイトを除く) をバイナリ形式でNSSファイルの末尾に出力
		print NSS substr($seq_sketch_index, 8);
		
		# NSSファイルを閉じる
		close(NSS);
		
		# スレッドを終了
		return(1);
	};
	## ここまでデータストリームスレッドの処理 ##
	
	# 変数を宣言
	my $read_id = 0;
	my $seq_read_index = "";
	
	# Nanoha Sequence Read (NSR) ファイルを作成
	open(NSR, ">", "$opt{o}.nsr") or &exception::error("failed to make file: $opt{o}.nsr");
	
	# NSRファイルをバイナリモードにする
	binmode(NSR);
	
	# ファイルヘッダーをバイナリ形式でNSRファイルに出力
	print NSR $file_header;
	
	# NSRファイルのファイルポインタの位置をシーケンスインデックスに登録
	vec($seq_read_index, 0, 64) = tell(NSR);
	
	# 入力の改行コードを一時的に変更 (-w指定時)
	local $/ = "\r\n" if $opt{"w"};
	
	# FASTQファイルを読み込みながら処理
	print STDERR "Sketching out sequence reads...";
	while (my $line = <>) {
		# 改行コードを除去
		chomp($line);
		
		# 読み込み行がシーケンス行以外の場合は以下の処理をスキップ
		next if $. % 4 != 2;
		
		# リードシーケンスを大文字に変換しリード両端を指定した長さだけトリミングしてシーケンスが消失した場合は以下の処理をスキップ
		my $read_seq = substr(uc($line), $opt{"5"}, length($line) - $opt{"5"} - $opt{"3"}) or next;
		
		# ACGT以外の塩基を除去
		$read_seq =~ s/[^ACGT]//g;
		
		# リード長を取得
		my $read_len = length($read_seq);
		
		# リード長がk-mer未満の場合は以下の処理をスキップ
		next if $read_len < $opt{"k"};
		
		# リードIDを更新
		$read_id++;
		
		# 実行中のスレッド数が指定値+1と一致している場合はリードIDとリードシーケンスをバイナリ形式で入力キューに追加
		$input->enqueue(pack("LA*", $read_id, $read_seq)) if threads->list(threads::running) == $opt{"p"} + 1;
		
		# リードシーケンスを2進数に変換
		$read_seq =~ s/A/00/g;
		$read_seq =~ s/C/10/g;
		$read_seq =~ s/G/01/g;
		$read_seq =~ s/T/11/g;
		
		# リード長とリードシーケンスをバイナリ形式でNSRファイルに出力
		print NSR pack("Nb*", $read_len, $read_seq);
		
		# NSRファイルのファイルポインタの位置をシーケンスインデックスに登録
		vec($seq_read_index, $read_id, 64) = tell(NSR);
		
		# リードIDが上限値に達した場合、読み込みを終了
		last if $read_id == $opt{"u"};
	}
	
	# シーケンスリードインデックスをバイナリ形式でNSRファイルの末尾に出力
	print NSR $seq_read_index;
	
	# NSRファイルを閉じる
	close(NSR);
	
	# 入力キューを終了
	$input->end;
	
	# 各ワーカースレッドが終了するまで待機
	$worker_threads[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 出力キューを終了
	$output->end;
	
	# データストリームスレッドが終了するまで待機
	$worker_threads[$opt{"p"}]->join or &exception::error("data stream thread abnormally exited");
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	return(1);
}

# weighted minimizerを生成 sketch::generate_weighted_minimizers(基数変換したk-merリストリファレンス, minimizer生成個数)
sub generate_weighted_minimizers {
	# 引数を取得
	my ($converted_kmers, $num_minimizers) = @_;
	
	# XSコードが有効の場合はXSの関数でminimizerを生成して返す
	return(&main::generate_weighted_minimizers($converted_kmers, $num_minimizers)) if $ENV{"NANOHA_XS"};
	
	# 変数を宣言
	my %hash_values = ();
	my $mt = undef;
	my $count = 0;
	
	# 基数変換したk-merをシード値に用いて指定した個数のハッシュ値を生成
	$mt = Math::Random::MT->new($_->[0]) and $count = $_->[1] and push(@{$hash_values{$_->[0]}}, map {1 - $mt->rand ** (1 / $count)} 1..$num_minimizers) foreach List::Util::pairs(@{$converted_kmers});
	
	# minimizerを生成
	my @minimizers = map {List::Util::reduce {$hash_values{$a}->[$_] < $hash_values{$b}->[$_] ? $a : $b} keys(%hash_values)} 0..$num_minimizers - 1;
	
	# minimizerを返す
	return(\@minimizers);
}

## ここからbuildコマンドのパッケージ ##
package build;

# コマンドとオプションを定義
sub define {
	$note = "Build sequence similarity graph under specified conditions.";
	$usage = "<prefix>";
	$option{"L PATH "} = "Path to a positive list file of k-mer minimizers";
	$option{"l PATH "} = "Path to a negative list file of k-mer minimizers";
	$option{"p INT "} = "Number of parallel worker threads <1-> [1]";
	$option{"m INT "} = "Maximum number of k-mer minimizers to be used from each sequence read <1-255> [255]";
	$option{"n INT "} = "Cutoff number of matched k-mer minimizers <1-255> [1]";
	$option{"v"} = "\tUse SIMD-vectorized code for LLCS calculation under XS code enabled";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("-L and -l incompatible: -L $opt{L}, -l $opt{l}") if defined($opt{"L"}) and defined($opt{"l"});
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT 1-255: -m $opt{m}") if $opt{"m"} !~ /^\d+$/ or $opt{"m"} < 1 or $opt{"m"} > 255;
	&exception::error("specify INT 1-255: -n $opt{n}") if $opt{"n"} !~ /^\d+$/ or $opt{"n"} < 1 or $opt{"n"} > 255;
	
	# 使用コードを表示
	print STDERR "Use ", $ENV{"NANOHA_XS"} ? "XS" : "Perl", " code\n";
	
	# 使用パラメータを表示
	print STDERR "Use positive list of k-mer minimizers: $opt{L}\n" if defined($opt{"L"});
	print STDERR "Use negative list of k-mer minimizers: $opt{l}\n" if defined($opt{"l"});
	print STDERR "Maximum number of k-mer minimizers to be used: $opt{m}\n";
	print STDERR "Cutoff number of matched k-mer minimizers: $opt{n}\n";
	
	# 入力ファイルのプレフィックスを取得
	my $input_prefix = shift(@ARGV);
	
	# 入力ファイルのプレフィックスを確認
	&exception::error("input file prefix not specified") unless defined($input_prefix);
	
	# 入力ファイルを確認
	&common::check_files(["$input_prefix.nsr", "$input_prefix.nss"]);
	
	# minimizersリストファイルを取得
	my $minimizers_list_file = defined($opt{"L"}) ? $opt{"L"} : defined($opt{"l"}) ? $opt{"l"} : undef;
	
	# 変数を宣言
	my @minimizers_list = ();
	
	# minimizersリストを定義 (minimizersリストファイル指定時)
	if (defined($minimizers_list_file)) {
		# minimizersリストファイルを確認
		&common::check_files([$minimizers_list_file]);
		
		# minimizersリストファイルを開く
		open(MINIMIZERS_LIST, "<", $minimizers_list_file) or &exception::error("failed to open file: $minimizers_list_file");
		
		# minimizersリストファイルを読み込みながら処理
		print STDERR "Loading minimizers list...";
		while (my $line = <MINIMIZERS_LIST>) {
			# 改行コードを除去
			chomp($line);
			
			# タブ文字で区切る
			my ($minimizer, $number) = split(/\t/, $line);
			
			# minimizerの配列を確認
			&exception::error("invalid minimizer: $minimizer") unless $minimizer =~ /^[ACGTacgt]{7,15}$/;
			
			# minimizerの配列を基数変換してminimizersリストに追加
			$minimizers_list[$number - 1]->{$convert->to_dec("C" . $minimizer)} = undef;
		}
		print STDERR "completed\n";
		
		# minimizersリストファイルを閉じる
		close(MINIMIZERS_LIST);
	}
	
	# Nanoha Sequence Read (NSR) ファイルを開く
	open(NSR, "<", "$input_prefix.nsr") or &exception::error("failed to open file: $input_prefix.nsr");
	
	# NSRファイルをバイナリモードにする
	binmode(NSR);
	
	# NSRファイルのファイルヘッダーを読み込む
	read(NSR, my $nsr_file_header, 8);
	
	# NSRファイルを閉じる
	close(NSR);
	
	# Nanoha Sequence Sketch (NSS) ファイルを開く
	open(NSS, "<", "$input_prefix.nss") or &exception::error("failed to open file: $input_prefix.nss");
	
	# NSSファイルをバイナリモードにする
	binmode(NSS);
	
	# NSSファイルのファイルヘッダーを読み込む
	read(NSS, my $nss_file_header, 8);
	
	# NSSファイルからパラメータ情報を読み込む
	read(NSS, my $parameter_info, 2);
	
	# ファイルバージョンの互換性を確認
	&exception::error("file version inconsistent with software version: $input_prefix.nss") if vec($parameter_info, 1, 4) >> 1 != [split(/\./, $version)]->[1] - 1;
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsr, $input_prefix.nss") unless $nsr_file_header eq $nss_file_header;
	
	# ストランド特異性を取得
	my $strand_speicific = vec($parameter_info, 4, 1);
	
	# minimizer数を定義
	my $num_minimizers = List::Util::min(vec($parameter_info, 1, 8), $opt{"m"});
	
	# リード数を取得
	my $num_reads = ((-s "$input_prefix.nss") - 10) / (vec($parameter_info, 1, 8) + 2) / 4;
	
	# シーケンススケッチインデックスを取得
	seek(NSS, -($num_reads + 1) * 8, 2);
	read(NSS, my $seq_sketch_index, ($num_reads + 1) * 8);
	
	# 変数を宣言
	my %minimizer_groups = ();
	
	# 各リードについて処理
	print STDERR "Loading sequence sketches...";
	for (my $read_id = 1;$read_id <= $num_reads;$read_id++) {
		# NSSファイルのファイルポインタをシーケンススケッチインデックスの位置にセット
		seek(NSS, vec($seq_sketch_index, $read_id, 64), 0);
		
		# NSSファイルからminimizer数×4バイト読み込む
		read(NSS, my $minimizers, $num_minimizers * 4);
		
		# ポジティブリストに該当しないminimizerを有する場合は除外
		next if defined($opt{"L"}) and List::MoreUtils::notall {exists($minimizers_list[$_]->{vec($minimizers, $_, 32)})} 0..$num_minimizers - 1;
		
		# ネガティブリストに該当するminimizerを有する場合は除外
		next if defined($opt{"l"}) and List::MoreUtils::any {exists($minimizers_list[$_]->{vec($minimizers, $_, 32)})} 0..$num_minimizers - 1;
		
		# minimizerグループにリードIDをバイナリ形式で追加
		$minimizer_groups{$minimizers} .= pack("L", $read_id);
	}
	print STDERR "completed\n";
	
	# NSSファイルを閉じる
	close(NSS);
	
	# 変数を宣言
	my @minimizer_buckets = map {{}} 1..$num_minimizers;
	
	# 各minimizersについて個々のminimizerごとにバケットソート
	print STDERR "Sorting k-mer minimizers...";
	List::MoreUtils::pairwise {$a->{$b} .= $_} @minimizer_buckets, @{[unpack("(a4)" . $num_minimizers, $_)]} foreach keys(%minimizer_groups);
	print STDERR "completed\n";
	
	### minimizerの一致数を算出 ###
	# Nanoha Sequence Graph (NSG) ファイルを作成
	open(NSG, ">", "$input_prefix.nsg") or &exception::error("failed to make file: $input_prefix.nsg");
	
	# NSGファイルをバイナリモードにする
	binmode(NSG);
	
	# ファイルヘッダーをバイナリ形式でNSGファイルに出力
	print NSG $nss_file_header;
	
	# 変数を宣言
	my @seq_graph_index = (0) x ($num_reads + 1);
	
	# 各minimizersについて処理
	print STDERR "Counting matched k-mer minimizers...";
	foreach my $minimizers (keys(%minimizer_groups)) {
		# minimizer一致数ごとにminimizersをバケットソート
		my $sorted_minimizers = &sort_minimizers(\@minimizer_buckets, $minimizers, $opt{"n"});
		
		# minimizer一致数ごとにリードIDを取得
		my @read_ids = map {join("", @minimizer_groups{@{$_}})} @{$sorted_minimizers};
		
		# minimizerが完全一致する各リードについてNSGファイルのファイルポインタの位置をシーケンスグラフインデックスに登録
		$seq_graph_index[$_] = tell(NSG) foreach unpack("L*", $read_ids[-1]);
		
		# minimizer一致数ごとのリード数とリードIDリストをバイナリ形式でNSGファイルに出力
		print NSG pack("L*", map {length($_) / 4} @read_ids), join("", @read_ids);
	}
	print STDERR "completed\n";
	
	# NSGファイルを閉じる
	close(NSG);
	
	# 変数を削除
	print STDERR "Releasing memories...";
	undef(%minimizer_groups);
	undef(@minimizer_buckets);
	print STDERR "completed\n";
	
	### LLCSに基づくリード間のJaccard類似度を算出 ###
	# 追加ヌル文字数を算出
	my $num_additional_nulls = $ENV{"NANOHA_XS"} && $opt{"v"} ? (0x10 << &main::get_word_size_index) - 1 : 0x0F;
	
	# 変数を宣言
	my @worker_threads = ();
	my $num_error_threads = 0;
	
	# 入出力キューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 入出力キューの要素数上限を定義
	$input->limit = $opt{"p"};
	$output->limit = $opt{"p"};
	
	# 指定されたワーカースレッド数で並列処理
	print STDERR "Calculating Jaccard similarities based on LLCS...";
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_threads[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# NSGファイルを開く
			open(NSG, "<", "$input_prefix.nsg") or &exception::error("failed to open file: $input_prefix.nsg");
			
			# NSGファイルをバイナリモードにする
			binmode(NSG);
			
			# NSRファイルを開く
			open(NSR, "<", "$input_prefix.nsr") or &exception::error("failed to open file: $input_prefix.nsr");
			
			# NSRファイルをバイナリモードにする
			binmode(NSR);
			
			# シーケンスリードインデックスを取得
			seek(NSR, -($num_reads + 1) * 8, 2);
			read(NSR, my $seq_read_index, ($num_reads + 1) * 8);
			
			# 入力キューを読み込みながら処理
			while (defined(my $query_read_id = $input->dequeue)) {
				# 変数を宣言
				my @read_seqs = ();
				my @read_lens = ();
				
				# NSGファイルのファイルポインタをシーケンスグラフインデックスの位置にセット
				seek(NSG, $seq_graph_index[$query_read_id], 0);
				
				# NSGファイルから(minimizer数-minimizer一致数閾値+1)×4バイト読み込む
				read(NSG, my $num_reads, ($num_minimizers - $opt{"n"} + 1) * 4);
				
				# minimizer一致数ごとのリード数を取得
				my @num_reads = unpack("L*", $num_reads);
				
				# NSGファイルからリード数×4バイト読み込む
				read(NSG, my $read_ids, List::Util::sum(@num_reads) * 4);
				
				# 各リードについて処理
				foreach my $read_id ($query_read_id, unpack("L*", $read_ids)) {
					# NSRファイルのファイルポインタをシーケンスリードインデックスの位置にセット
					seek(NSR, vec($seq_read_index, $read_id - 1, 64), 0);
					
					# NSRファイルから4バイト読み込みリード長リストに追加
					read(NSR, my $read_len, 4);
					
					# NSRファイルからリードシーケンスをバイナリ形式で読み込む
					read(NSR, my $read_seq, vec($seq_read_index, $read_id, 64) - tell(NSR));
					
					# リードシーケンスリストにリードシーケンスをバイナリ形式で追加
					push(@read_seqs, $read_seq . chr(0) x $num_additional_nulls);
					
					# リード長リストにリード長を追加
					push(@read_lens, vec($read_len, 0, 32));
				}
				
				# LLCSに基づくJaccard類似度を算出
				my $LLCS_Jaccard_similarities = &calc_LLCS_Jaccard_similarities(\@read_seqs, \@read_lens, $strand_speicific, $opt{"v"});
				
				# 各リードについてminimizer一致数をリスト化
				my @num_matched_minimizers = map {($_ + $opt{"n"}) x $num_reads[$_]} 0..$num_minimizers - $opt{"n"};
				
				# LLCSに基づくJaccard類似度で修正した重みをクエリーリードID及びリード数とともにバイナリ形式で出力キューに追加
				$output->enqueue(pack("LLF*", $query_read_id, scalar(@num_matched_minimizers), List::MoreUtils::pairwise {$a * $b / $num_minimizers} @{$LLCS_Jaccard_similarities}, @num_matched_minimizers));
			}
			
			# NSRファイルを閉じる
			close(NSR);
			
			# NSGファイルを閉じる
			close(NSG);
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# データストリームスレッドを作成
	## ここからデータストリームスレッドの処理 ##
	$worker_threads[$opt{"p"}] = threads::async {
		# このスレッドのみを終了可能に変更
		threads->set_thread_exit_only(1);
		
		# NSGファイルを追加書き込みモードで開く
		open(NSG, ">>", "$input_prefix.nsg") or &exception::error("failed to open file: $input_prefix.nsg");
		
		# NSGファイルをバイナリモードにする
		binmode(NSG);
		
		# 出力キューを読み込みながら処理
		while (defined(my $serialized_data = $output->dequeue)) {
			# クエリーリードID及びリード数を取得
			my ($query_read_id, $num_reads) = unpack("LL", $serialized_data);
			
			# リレーするシーケンスグラフインデックスを算出
			my $relayed_seq_graph_index = $seq_graph_index[$query_read_id] + ($num_minimizers - $opt{"n"} + 1) * 4;
			
			# NSGファイルのファイルポインタの位置をシーケンスグラフインデックスに登録して更新
			$seq_graph_index[$query_read_id] = tell(NSG);
			
			# リレーするシーケンスグラフインデックス及び重みをバイナリ形式でNSGファイルに出力
			print NSG pack("Q", $relayed_seq_graph_index), $serialized_data;
		}
		
		# 更新後のシーケンスグラフインデックスをバイナリ形式でNSGファイルに出力
		print NSG pack("Q*", @seq_graph_index[1..$#seq_graph_index]);
		
		# NSGファイルを閉じる
		close(NSG);
		
		# スレッドを終了
		return(1);
	};
	## ここまでデータストリームスレッドの処理 ##
	
	# シーケンスグラフインデックスが定義されている各リードについて実行中のスレッド数が指定値+1と一致している場合はリードを入力キューに追加
	threads->list(threads::running) == $opt{"p"} + 1 and $input->enqueue($_) foreach grep {$seq_graph_index[$_]} 1..$num_reads;
	
	# 入力キューを終了
	$input->end;
	
	# 各ワーカースレッドが終了するまで待機
	$worker_threads[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 出力キューを終了
	$output->end;
	
	# データストリームスレッドが終了するまで待機
	$worker_threads[$opt{"p"}]->join or &exception::error("data stream thread abnormally exited");
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	# 抽出したリード数を取得
	my $num_extracted_reads = List::MoreUtils::true {$_} @seq_graph_index;
	
	# 抽出したリード数を表示
	print STDERR $num_extracted_reads ? $num_extracted_reads : "no", " sequence read", $num_extracted_reads == 1 ? "" : "s", " extracted\n";
	
	return(1);
}

# minimizer一致数ごとにminimizersをバケットソート build::sort_minimizers(minimizerバケットリストリファレンス, minimizers, minimizer一致数閾値)
sub sort_minimizers {
	# 引数を取得
	my ($minimizer_buckets, $minimizers, $cutoff_num_matched_minimizers) = @_;
	
	# XSコードが有効の場合はXSの関数でminimizer一致数ごとにminimizersをバケットソートして返す
	return(&main::sort_minimizers($minimizer_buckets, $minimizers, $cutoff_num_matched_minimizers)) if $ENV{"NANOHA_XS"};
	
	# minimizer一致数ごとにminimizersをバケットソート
	my @sorted_minimizers = List::MoreUtils::occurrences map {unpack("(a" . length($minimizers) . ")*", $_)}
	grep {defined($_)} List::MoreUtils::pairwise {$a->{$b}} @{$minimizer_buckets}, @{[unpack("(a4)" . scalar(@{$minimizer_buckets}), $minimizers)]};
	
	# minimizer一致数が指定した個数未満のminimizersを削除
	splice(@sorted_minimizers, 0, $cutoff_num_matched_minimizers);
	
	# 未定義値を空リストリファレンスに変換して返す
	return([map {defined($_) ? $_ : []} @sorted_minimizers]);
}

# LLCSに基づくJaccard類似度を算出 build::calc_LLCS_Jaccard_similarities(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性, ベクトル演算コード有効フラグ)
sub calc_LLCS_Jaccard_similarities {
	# 引数を取得
	my ($read_seqs, $read_lens, $strand_speicific, $vectorized) = @_;
	
	# XSコードが有効の場合はXSの関数でLLCSに基づくJaccard類似度を算出して返す
	return($vectorized ? &main::calc_LLCS_Jaccard_similarities_by_simd($read_seqs, $read_lens, $strand_speicific) : &main::calc_LLCS_Jaccard_similarities($read_seqs, $read_lens, $strand_speicific)) if $ENV{"NANOHA_XS"};
	
	# クエリーリードシーケンスを取得
	my $query_seq = [unpack("Q*", $read_seqs->[0])];
	
	# クエリーリード長を取得
	my $query_len = $read_lens->[0];
	
	# クエリーブロック数を算出
	my $num_query_blocks = ($query_len >> 6) + (($query_len & 0x3F) > 0);
	
	# 変数を宣言
	my @LLCS_Jaccard_similarities = ();
	my @query_matrix = ([], [], [], []);
	my @query_base = ();
	my $mask = 0xFFFFFFFFFFFFFFFF;
	
	# クエリー行列を作成
	for (my $i = 0;$i < @query_matrix;$i++) {
		for (my $j = 0;$j < $num_query_blocks;$j++) {
			$query_base[0] = $query_seq->[$j * 2] ^ $mask;
			$query_base[0] &= $query_base[0] >> 1;
			$query_base[1] = $query_seq->[$j * 2 + 1] ^ $mask;
			$query_base[1] &= $query_base[1] << 1;
			$query_matrix[$i]->[$j] = 0;
			$query_matrix[$i]->[$j] |= $query_base[0] >> $_ & 0x0000000000000001 << $_ | $query_base[1] << $_ & 0x8000000000000000 >> $_ foreach 0..31;
		}
		$mask -= 0x5555555555555555;
	}
	
	# クエリー行列の端をマスク
	$query_matrix[0]->[$num_query_blocks - 1] &= 0xFFFFFFFFFFFFFFFF >> (-$query_len & 0x3F);
	
	# 残りの各リードについて処理
	for (my $n = 1;$n < @{$read_seqs};$n++) {
		# リードシーケンスを取得
		my $read_seq = [unpack("Q*", $read_seqs->[$n])];
		
		# リード長を取得
		my $read_len = $read_lens->[$n];
		
		# ブロック数を算出
		my $num_blocks = ($read_len >> 5) + (($read_len & 0x1F) > 0);
		
		# 変数を宣言
		my @llcs = ();
		my @v = ();
		my $u = 0;
		my $w = 0;
		
		# 順鎖及び逆鎖について処理
		for (my $strand = 0;$strand < 2;$strand++) {
			# 計算要素を初期化
			@v = (0xFFFFFFFFFFFFFFFF) x $num_query_blocks;
			
			# LLCSを算出
			for (my $i = $strand * ($num_blocks - 1);$i >= 0 and $i < $num_blocks;$i += 1 - 2 * $strand) {
				my $end = $i + !!($read_len & 0x1F) < $num_blocks ? 64 : $read_len << 1 & 0x3E;
				for (my $j = $strand * ($end - 2);$j >= 0 and $j < $end;$j += 2 - 4 * $strand) {
					my $base = $read_seq->[$i] >> $j & 0x03;
					my $carry = 0;
					for (my $k = 0;$k < $num_query_blocks;$k++) {
						$u = $v[$k] & $query_matrix[$base]->[$k];
						$w = ($v[$k] & 0x7FFFFFFFFFFFFFFF) + ($u & 0x7FFFFFFFFFFFFFFF) + $carry;
						$w ^= ($v[$k] ^ $u) & 0x8000000000000000;
						$carry = $w < $v[$k];
						$v[$k] = $w | $v[$k] - $u;
					}
				}
			}
			$llcs[$strand] = List::Util::sum(
				map {($_ & 0x00000000FFFFFFFF) + ($_ >> 32 & 0x00000000FFFFFFFF)}
				map {($_ & 0x0000FFFF0000FFFF) + ($_ >> 16 & 0x0000FFFF0000FFFF)}
				map {($_ & 0x00FF00FF00FF00FF) + ($_ >> 8 & 0x00FF00FF00FF00FF)}
				map {($_ & 0x0F0F0F0F0F0F0F0F) + ($_ >> 4 & 0x0F0F0F0F0F0F0F0F)}
				map {($_ & 0x3333333333333333) + ($_ >> 2 & 0x3333333333333333)}
				map {($_ & 0x5555555555555555) + ($_ >> 1 & 0x5555555555555555)}
				map {~$_} @v
			);
			
			# ストランド特異性が有効な場合はループを終了
			last if $strand_speicific;
			
			# リードシーケンスを相補的な塩基に変換
			$_ = ~$_ foreach @{$read_seq};
		}
		# ストランド特異性が有効な場合は順鎖に対するLLCSに基づくJaccard類似度を算出して以下の処理をスキップ
		push(@LLCS_Jaccard_similarities, $llcs[0] / ($query_len + $read_len - $llcs[0])) and next if $strand_speicific;
		
		# LLCSに基づくJaccard類似度を算出 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		push(@LLCS_Jaccard_similarities, (1 - 2 * ($llcs[0] < $llcs[1])) * $llcs[$llcs[0] < $llcs[1]] / ($query_len + $read_len - $llcs[$llcs[0] < $llcs[1]]));
	}
	
	# LLCSに基づくJaccard類似度を返す
	return(\@LLCS_Jaccard_similarities);
}

## ここからassortコマンドのパッケージ ##
package assort;

# コマンドとオプションを定義
sub define {
	$note = "Assort sequence reads based on sequence similarity graph under specified conditions.";
	$usage = "<prefix>";
	$option{"p INT "} = "Number of parallel worker threads <1-> [1]";
	$option{"t INT "} = "Number of trials per worker thread for clustering sequence reads <1-> [1]";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT >= 1: -t $opt{t}") if $opt{"t"} !~ /^\d+$/ or $opt{"t"} < 1;
	
	# 使用コードを表示
	print STDERR "Use ", $ENV{"NANOHA_XS"} ? "XS" : "Perl", " code\n";
	
	# 使用パラメータを表示
	print STDERR "Number of trials per worker thread for clustering: $opt{t}\n";
	
	# 入力ファイルのプレフィックスを取得
	my $input_prefix = shift(@ARGV);
	
	# 入力ファイルのプレフィックスを確認
	&exception::error("input file prefix not specified") unless defined($input_prefix);
	
	# 入力ファイルを確認
	&common::check_files(["$input_prefix.nss", "$input_prefix.nsg"]);
	
	# Nanoha Sequence Sketch (NSS) ファイルを開く
	open(NSS, "<", "$input_prefix.nss") or &exception::error("failed to open file: $input_prefix.nss");
	
	# NSSファイルをバイナリモードにする
	binmode(NSS);
	
	# NSSファイルのファイルヘッダーを読み込む
	read(NSS, my $nss_file_header, 8);
	
	# NSSファイルからパラメータ情報を読み込む
	read(NSS, my $parameter_info, 2);
	
	# NSSファイルを閉じる
	close(NSS);
	
	# ファイルバージョンの互換性を確認
	&exception::error("file version inconsistent with software version: $input_prefix.nss") if vec($parameter_info, 1, 4) >> 1 != [split(/\./, $version)]->[1] - 1;
	
	# ストランド特異性を取得
	my $strand_speicific = vec($parameter_info, 4, 1);
	
	# リード数を取得
	my $num_reads = ((-s "$input_prefix.nss") - 10) / (vec($parameter_info, 1, 8) + 2) / 4;
	
	# Nanoha Sequence Graph (NSG) ファイルを開く
	open(NSG, "<", "$input_prefix.nsg") or &exception::error("failed to open file: $input_prefix.nsg");
	
	# NSGファイルをバイナリモードにする
	binmode(NSG);
	
	# NSGファイルのファイルヘッダーを読み込む
	read(NSG, my $nsg_file_header, 8);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsg, $input_prefix.nss") unless $nsg_file_header eq $nss_file_header;
	
	# シーケンスグラフインデックスを取得
	seek(NSG, -($num_reads + 1) * 8, 2);
	read(NSG, my $seq_graph_index, ($num_reads + 1) * 8);
	my @seq_graph_index = unpack("Q*", $seq_graph_index);
	
	# 変数を宣言
	my @max_clusters = ();
	my $already_cluster_assigned = "";
	
	print STDERR "Creating maximum sequence clusters...";
	# 各リードについて到達可能なリードを列挙し極大クラスターを決定
	for (my $read_id = 1;$read_id <= $num_reads;$read_id++) {
		# シーケンスグラフインデックスが0の場合は以下の処理をスキップ
		next unless $seq_graph_index[$read_id];
		
		# クラスターに割り当て済みの場合は以下の処理をスキップ
		next if vec($already_cluster_assigned, $read_id, 1);
		
		# 変数を宣言
		my @read_queue = ($read_id);
		
		# クラスターを作成
		my $max_cluster = pack("L", $read_id);
		
		# クラスター割り当てフラグをチェック
		vec($already_cluster_assigned, $read_id, 1) = 1;
		
		# リードキューからリードIDを取得して処理
		while (my $read_id = shift(@read_queue)) {
			# シーケンスグラフファイルのファイルポインタをシーケンスグラフインデックスの位置にセット
			seek(NSG, $seq_graph_index[$read_id], 0);
			
			# シーケンスグラフファイルから16バイト読み込む
			read(NSG, my $serialized_data, 16);
			
			# リレーするシーケンスグラフインデックス及びリード数を取得
			my ($relayed_seq_graph_index, $num_reads) = unpack("Qx4L", $serialized_data);
			
			# シーケンスグラフファイルのファイルポインタをリレーするシーケンスグラフインデックスの位置にセット
			seek(NSG, $relayed_seq_graph_index, 0);
			
			# シーケンスグラフファイルからリード数×4バイト読み込む
			read(NSG, $serialized_data, $num_reads * 4);
			
			# リードIDを取得
			my @read_ids = unpack("L*", $serialized_data);
			
			# クラスターに未割り当てのリードについてリードキューとクラスターに追加してクラスター割り当てフラグをチェック
			push(@read_queue, $_) and $max_cluster .= pack("L", $_) and vec($already_cluster_assigned, $_, 1) = 1 foreach grep {!vec($already_cluster_assigned, $_, 1)} @read_ids;
		}
		
		# 作成したクラスターをバイナリ形式で極大クラスターリストに追加
		push(@max_clusters, $max_cluster);
	}
	print STDERR "completed\n";
	
	# NSGファイルを閉じる
	close(NSG);
	
	### excess modularity density (Qx) に基づくリードのクラスタリング ###
	# 共有変数を宣言
	my $max_Qx : shared;
	my $best_cluster_assignment : shared;
	
	# Qx最大値の初期値を設定
	$max_Qx = -Inf;
	
	# 変数を宣言
	my @worker_threads = ();
	my $num_error_threads = 0;
	
	# 指定されたワーカースレッド数で並列処理
	print STDERR "Creating sequence clusters...";
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_threads[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 試行回数だけ処理
			for (my $i = 0;$i < $opt{"t"};$i++) {
				# クラスター割り当てリストを初期化
				my $cluster_assignment = pack("L*", 0..$num_reads);
				
				# excess modularity density (Qx) に基づいてリードをクラスタリング
				my $Qx = &create_sequence_clusters(\$cluster_assignment, \@max_clusters, $seq_graph_index, "$input_prefix.nsg");
				
				# QxにNaNが返ってきた場合はファイルオープンのエラーとする
				&exception::error("failed to open file: $input_prefix.nsg") unless $Qx == $Qx;
				
				# 共有変数をロック
				lock($max_Qx);
				
				# Qxが最大値を超えた場合はQx最大値及び最善のクラスター割り当てリストを更新
				($max_Qx, $best_cluster_assignment) = ($Qx, $cluster_assignment) if $Qx > $max_Qx;
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 各ワーカースレッドが終了するまで待機
	$worker_threads[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	### 判別分析法（大津の二値化）に基づく深度の推奨閾値を算出 ###
	# 変数を宣言
	my %clusters = ();
	my %small_clusters = ();
	my %large_clusters = ();
	my %cutoff_depth = ();
	my $num_small_clusters = 0;
	my $num_large_clusters = 0;
	my $sum_small_clusters = 0;
	my $sum_large_clusters = 0;
	
	# 各リードについて所属するクラスターに登録
	print STDERR "Determining cutoff cluster size...";
	$clusters{substr($best_cluster_assignment, $_ * 4, 4)} .= pack("N", $_) foreach grep {$seq_graph_index[$_]} 1..$num_reads;
	
	# クラスターサイズごとにクラスター数を算出
	$small_clusters{$_}++ foreach map {length($_) / 4} values(%clusters);
	
	# 閾値を変えながら分離度の指標値を算出
	$large_clusters{$_} = delete($small_clusters{$_}) and %small_clusters and
	($num_small_clusters, $num_large_clusters) = (List::Util::sum(values(%small_clusters)), List::Util::sum(values(%large_clusters))) and
	($sum_small_clusters, $sum_large_clusters) = (List::Util::sum(List::Util::pairmap {$a * $b} %small_clusters), List::Util::sum(List::Util::pairmap {$a * $b} %large_clusters)) and
	$cutoff_depth{$_} = $num_small_clusters * $num_large_clusters * ($sum_small_clusters / $num_small_clusters - $sum_large_clusters / $num_large_clusters) ** 2
	foreach sort {$b <=> $a} keys(%small_clusters);
	
	# 深度の推奨閾値を算出
	my $recommended_cutoff_depth = List::Util::reduce {$cutoff_depth{$a} > $cutoff_depth{$b} ? $a : $cutoff_depth{$a} == $cutoff_depth{$b} && $a < $b ? $a : $b} keys(%cutoff_depth);
	print STDERR "completed\n";
	
	### シーケンス方向の決定 ###
	# 変数を宣言
	my $seq_ori = pack("C", $strand_speicific);
	
	# Nanoha Sequence Cluster (NSC) ファイルを作成
	open(NSC, ">", "$input_prefix.nsc") or &exception::error("failed to make file: $input_prefix.nsc");
	
	# NSCファイルをバイナリモードにする
	binmode(NSC);
	
	# ファイルヘッダーと深度の推奨閾値をバイナリ形式でNSCファイルに出力
	print NSC $nss_file_header, pack("N", $recommended_cutoff_depth);
	
	# Nanoha Sequence Graph (NSG) ファイルを開く
	open(NSG, "<", "$input_prefix.nsg") or &exception::error("failed to open file: $input_prefix.nsg");
	
	# NSGファイルをバイナリモードにする
	binmode(NSG);
	
	# 各クラスターについて処理
	print STDERR "Determining sequence orientations...";
	foreach my $cluster (values(%clusters)) {
		# 所属リードのリードIDをバイナリ形式でNSCファイルに出力
		print NSC pack("N", length($cluster) / 4), $cluster;
		
		# 所属リードのリードIDを取得してリードキューに追加
		my @read_queue = unpack("N*", $cluster);
		
		# リードキューから先頭リードを除いた残りの所属リードのリードIDをメンバーに登録
		my %members = map {$_ => 1} splice(@read_queue, 1);
		
		# リードキューからリードIDを取得して処理
		while (my $read_id = shift(@read_queue)) {
			# シーケンスグラフファイルのファイルポインタをシーケンスグラフインデックスの位置にセット
			seek(NSG, $seq_graph_index[$read_id], 0);
			
			# シーケンスグラフファイルから16バイト読み込む
			read(NSG, my $serialized_data, 16);
			
			# リレーするシーケンスグラフインデックス及びリード数を取得
			my ($relayed_seq_graph_index, $num_reads) = unpack("Qx4L", $serialized_data);
			
			# シーケンスグラフファイルからリード数×8バイト読み込む
			read(NSG, $serialized_data, $num_reads * 8);
			
			# シーケンスグラフファイルのファイルポインタをリレーするシーケンスグラフインデックスの位置にセット
			seek(NSG, $relayed_seq_graph_index, 0);
			
			# シーケンスグラフファイルからリード数×4バイト読み込む
			read(NSG, $serialized_data, $num_reads * 4, length($serialized_data));
			
			# リードIDと重みを取得
			my @read_ids = unpack("F" . $num_reads . "L" . $num_reads, $serialized_data);
			my @weights = splice(@read_ids, 0, $num_reads);
			
			# リンクする各リードのうち、メンバーに登録されているものついてメンバーから削除するとともにリードキューに追加してシーケンス方向を決定
			delete($members{$read_ids[$_]}) and push(@read_queue, $read_ids[$_]) and vec($seq_ori, $read_ids[$_], 1) = vec($seq_ori, $read_id, 1) ^ ($weights[$_] < 0) foreach 0..$num_reads - 1;
		}
	}
	print STDERR "completed\n";
	
	# NSGファイルを閉じる
	close(NSG);
	
	# NSCファイルを閉じる
	close(NSC);
	
	# Nanoha Sequence Orientation (NSO) ファイルを作成
	open(NSO, ">", "$input_prefix.nso") or &exception::error("failed to make file: $input_prefix.nso");
	
	# NSOファイルをバイナリモードにする
	binmode(NSO);
	
	# ファイルヘッダー及びシーケンス方向をバイナリ形式でNSOファイルに出力
	print NSO $nss_file_header, $seq_ori;
	
	# NSOファイルを閉じる
	close(NSO);
	
	# クラスター数を表示
	print STDERR scalar(keys(%clusters)), " sequence cluster", %clusters > 2 ? "s" : "", " found\n";
	
	# Qx最大値を表示
	print STDERR "Excess modularity density: $max_Qx\n";
	
	# 深度の推奨閾値を表示
	print STDERR "Recommended cutoff depth: $recommended_cutoff_depth\n" if defined($recommended_cutoff_depth);
	
	return(1);
}

# 密度を算出 assort::calc_density(内部エッジの重みの総和, ノード数)
sub calc_density {
	# 引数を取得
	my ($sum_internal_edges, $num_nodes) = @_;
	
	# 密度を算出
	my $density = $num_nodes > 1 ? $sum_internal_edges / $num_nodes / ($num_nodes - 1) : 0;
	
	# 密度を返す
	return($density);
}

# 部分的なQxを算出 assort::calc_partial_Qx(内部エッジの重みの総和, ノード数, クラスター次数, グラフ次数, グラフ密度)
sub calc_partial_Qx {
	# 引数を取得
	my ($sum_internal_edges, $num_nodes, $cluster_degree, $graph_degree, $graph_density) = @_;
	
	# 密度指数を算出
	my $density_index = &calc_density($sum_internal_edges, $num_nodes) - $graph_density;
	
	# 部分的なQxを算出
	my $partial_Qx = $sum_internal_edges * $density_index - ($cluster_degree * $density_index) ** 2 / $graph_degree;
	
	# 部分的なQxを返す
	return($partial_Qx);
}

# excess modularity density (Qx) に基づくリードのクラスタリング assort::create_sequence_clusters(クラスター割り当てリストリファレンス, 極大クラスターリストリファレンス, シーケンスグラフインデックス, nanohaシーケンスグラフファイルパス)
sub create_sequence_clusters {
	# 引数を取得
	my ($cluster_assignment, $max_clusters, $seq_graph_index, $nanoha_sequence_graph) = @_;
	
	# XSコードが有効の場合はXSの関数でexcess modularity density (Qx) に基づくリードのクラスタリングを実行して返す
	return(&main::create_sequence_clusters(${$cluster_assignment}, $max_clusters, $seq_graph_index, $nanoha_sequence_graph)) if $ENV{"NANOHA_XS"};
	
	# リード数を取得
	my $num_reads = length(${$cluster_assignment}) / 4 - 1;
	
	# シーケンスグラフインデックスを取得
	my @seq_graph_index = unpack("Q*", $seq_graph_index);
	
	# 極大クラスターを取得
	$max_clusters = [map {[unpack("L*", $_)]} @{$max_clusters}];
	
	# NSGファイルを開く
	open(NSG, "<", $nanoha_sequence_graph) or &exception::error("failed to open file: $nanoha_sequence_graph");
	
	# NSGファイルをバイナリモードにする
	binmode(NSG);
	
	# 変数を宣言
	my @node_degree = (0) x ($num_reads + 1);
	my @sum_internal_edges = (0) x ($num_reads + 1);
	my @num_nodes = (1) x ($num_reads + 1);
	my @cluster_degree = (0) x ($num_reads + 1);
	my @node_queue = ();
	my $cluster_updated = "";
	my $last_Qx = -Inf;
	
	# 変数を初期化
	for (my $read_id = 1;$read_id <= $num_reads;$read_id++) {
		# シーケンスグラフインデックスが0の場合は以下の処理をスキップ
		next unless $seq_graph_index[$read_id];
		
		# シーケンスグラフファイルのファイルポインタをシーケンスグラフインデックスの位置にセット
		seek(NSG, $seq_graph_index[$read_id], 0);
		
		# シーケンスグラフファイルから16バイト読み込む
		read(NSG, my $serialized_data, 16);
		
		# リード数を取得
		my ($num_reads) = unpack("x12L", $serialized_data);
		
		# シーケンスグラフファイルからリード数×8バイト読み込む
		read(NSG, $serialized_data, $num_reads * 8);
		
		# 重みの和を取得 (自己リンクの1を減算)
		$node_degree[$read_id] = List::Util::sum(map {abs($_)} unpack("F*", $serialized_data)) - 1;
		$cluster_degree[$read_id] = $node_degree[$read_id];
		
		# ノードキューを初期化
		push(@node_queue, $read_id);
	}
	
	# シーケンスグラフインデックスが0でないリード数を取得
	$num_reads = @node_queue;
	
	# クラスター割り当てリストを取得
	my @cluster_assignment = unpack("L*", ${$cluster_assignment});
	
	# グラフ次数を算出
	my $graph_degree = List::Util::sum(@node_degree);
	
	# グラフ密度を算出
	my $graph_density = &calc_density($graph_degree, $num_reads);
	
	# 初期Qxを算出
	my $Qx = List::Util::sum(map {&calc_partial_Qx(0, 1, $_, $graph_degree, $graph_density)} @cluster_degree);
	
	# メルセンヌ・ツイスターによる擬似乱数生成器を作成
	my $mt = Math::Random::MT->new;
	
	# Qxが向上した場合は処理を継続
	while ($Qx > $last_Qx) {
		# 現在のQxを保存
		$last_Qx = $Qx;
		
		# クラスター更新フラグをリセット
		$cluster_updated = "";
		
		# 各ノードについて処理
		for (my $i = $num_reads - 1;$i >= 0;$i--) {
			# ランダムにノードを選択
			my $j = int($mt->rand($i + 1));
			my $node = $node_queue[$j];
			$node_queue[$j] = $node_queue[$i];
			$node_queue[$i] = $node;
			
			# シーケンスグラフファイルのファイルポインタをシーケンスグラフインデックスの位置にセット
			seek(NSG, $seq_graph_index[$node], 0);
			
			# シーケンスグラフファイルから16バイト読み込む
			read(NSG, my $serialized_data, 16);
			
			# リレーするシーケンスグラフインデックス及びリード数を取得
			my ($relayed_seq_graph_index, $num_reads) = unpack("Qx4L", $serialized_data);
			
			# シーケンスグラフファイルからリード数×8バイト読み込む
			read(NSG, $serialized_data, $num_reads * 8);
			
			# シーケンスグラフファイルのファイルポインタをリレーするシーケンスグラフインデックスの位置にセット
			seek(NSG, $relayed_seq_graph_index, 0);
			
			# シーケンスグラフファイルからリード数×4バイト読み込む
			read(NSG, $serialized_data, $num_reads * 4, length($serialized_data));
			
			# リードIDと重みを取得
			my @read_ids = unpack("F" . $num_reads . "L" . $num_reads, $serialized_data);
			my @weights = map {abs($_)} splice(@read_ids, 0, $num_reads);
			
			# 所属クラスターを取得
			my $cluster = $cluster_assignment[$node];
			
			# 変数を宣言 (所属クラスターは自己リンクの2をあらかじめ減算)
			my %sum_linked_edges = ($cluster => -2, 0 => 0);
			
			# 隣接クラスターに対するエッジの重みの総和を算出
			$sum_linked_edges{$cluster_assignment[$read_ids[$_]]} += 2 * $weights[$_] foreach 0..$num_reads - 1;
			
			# 所属クラスターから転出する場合のQx変化量を算出
			my $basal_Qx_change = &calc_partial_Qx($sum_internal_edges[$cluster] - $sum_linked_edges{$cluster}, $num_nodes[$cluster] - 1, $cluster_degree[$cluster] - $node_degree[$node], $graph_degree, $graph_density) - &calc_partial_Qx($sum_internal_edges[$cluster], $num_nodes[$cluster], $cluster_degree[$cluster], $graph_degree, $graph_density);
			
			# 隣接クラスターに転入 (あるいは所属クラスターから独立) する場合のQx変化量を算出
			my %Qx_change = map {$_ => &calc_partial_Qx($sum_internal_edges[$_] + $sum_linked_edges{$_}, $num_nodes[$_] + 1, $cluster_degree[$_] + $node_degree[$node], $graph_degree, $graph_density) - &calc_partial_Qx($sum_internal_edges[$_], $num_nodes[$_], $cluster_degree[$_], $graph_degree, $graph_density)} keys(%sum_linked_edges);
			
			# 所属クラスターを変えない場合のQx変化量を定義
			$Qx_change{$cluster} = -$basal_Qx_change;
			
			# Qx変化量が最大となるクラスターを選出
			my $target_cluster = List::Util::reduce {$Qx_change{$a} > $Qx_change{$b} ? $a : $Qx_change{$a} == $Qx_change{$b} && $a > $b ? $a : $b} keys(%Qx_change);
			
			# 所属クラスターから独立する場合は新規クラスターを登録
			$target_cluster = List::MoreUtils::firstidx {!$_} @num_nodes and ($sum_linked_edges{$target_cluster}, $Qx_change{$target_cluster}) = (0, $Qx_change{0}) unless $target_cluster;
			
			# クラスター内エッジの重みの総和を更新
			$sum_internal_edges[$target_cluster] += $sum_linked_edges{$target_cluster};
			$sum_internal_edges[$cluster] -= $sum_linked_edges{$cluster};
			
			# ノード数を更新
			$num_nodes[$target_cluster]++;
			$num_nodes[$cluster]--;
			
			# クラスター次数を更新
			$cluster_degree[$target_cluster] += $node_degree[$node];
			$cluster_degree[$cluster] -= $node_degree[$node];
			
			# 所属クラスターを更新
			$cluster_assignment[$node] = $target_cluster;
			
			# Qxを更新
			$Qx += $Qx_change{$target_cluster} + $basal_Qx_change;
			
			# クラスター更新フラグをチェック
			vec($cluster_updated, $target_cluster, 1) = $target_cluster != $cluster;
			vec($cluster_updated, $cluster, 1) = $target_cluster != $cluster;
		}
		
		# 各極大クラスターについて属するノードの所属クラスターのいずれかが更新されたか否かで分類
		(my $fixed_clusters, $max_clusters) = List::MoreUtils::part {List::MoreUtils::any {vec($cluster_updated, $cluster_assignment[$_], 1)} @{$_}} @{$max_clusters};
		
		# 所属クラスターが固定されたノードを取得
		my %fixed_nodes = map {$_ => undef} map {@{$_}} @{$fixed_clusters};
		
		# 所属クラスターが固定された各ノードについてノードキューの位置を取得
		$fixed_nodes{$node_queue[$_]} = $_ foreach grep {exists($fixed_nodes{$node_queue[$_]})} 0..$num_reads - 1;
		
		# 固定されたノード数だけリード数を減算
		$num_reads -= keys(%fixed_nodes);
		
		# 変数を宣言
		my $fixed_node = undef;
		
		# 所属クラスターが固定された各ノードについてノードキューの末尾に配置
		$fixed_node = [each(%fixed_nodes)] and @node_queue[$fixed_node->[1], $_] = @node_queue[$_, $fixed_node->[1]] and
		exists($fixed_nodes{$node_queue[$fixed_node->[1]]}) and $fixed_nodes{$node_queue[$fixed_node->[1]]} = $fixed_node->[1]
		foreach $num_reads..$num_reads + keys(%fixed_nodes) - 1;
	}
	
	# NSGファイルを閉じる
	close(NSG);
	
	# クラスター割り当てリストをバイナリ形式に変換
	${$cluster_assignment} = pack("L*", @cluster_assignment);
	
	# Qxを返す
	return($Qx / $graph_degree);
}

## ここからunifyコマンドのパッケージ ##
package unify;

# コマンドとオプションを定義
sub define {
	$note = "Unify sequence reads in the same cluster under specified conditions.";
	$usage = "<prefix> [>out.fa]";
	$option{"c STR "} = "Alignment method to generate consensus sequence reads <local|global|semi-global>";
	$option{"d INT "} = "Cutoff depth (cluster size) to eliminate small clusters <1->";
	$option{"p INT "} = "Number of parallel worker threads <1-> [1]";
	$option{"g INT "} = "Gap opening penalty <1-> [3]";
	$option{"e INT "} = "Gap extension penalty <1-> [1]";
	$option{"m INT "} = "Match award <1-> [5]";
	$option{"n INT "} = "Mismatch penalty <0-> [4]";
	$option{"a INT "} = "Maximum amount of sequence reads to be aligned <1-" . main::max_depth . "> [" . main::max_depth . "]";
	$option{"q FLOAT "} = "Cutoff false discovery rate to eliminate strand-biased clusters <0-1> [0.001]";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("unknown algorithm specified: -c $opt{c}") if defined($opt{"c"}) and !grep {/^$opt{"c"}$/i} ("local", "global", "semi-global");
	&exception::error("generating consensus sequence disabled (NANOHA_XS = 0 or NANOHA_SPOA = 0): -c $opt{c}") if defined($opt{"c"}) and (!$ENV{"NANOHA_XS"} or !$ENV{"NANOHA_SPOA"});
	&exception::error("specify INT >= 1: -d $opt{d}") if defined($opt{"d"}) and ($opt{"d"} !~ /^\d+$/ or $opt{"d"} < 1);
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT >= 1: -g $opt{g}") if $opt{"g"} !~ /^\d+$/ or $opt{"g"} < 1;
	&exception::error("specify INT <= 1: -e $opt{e}") if $opt{"e"} !~ /^\d+$/ or $opt{"e"} < 1;
	&exception::error("specify INT >= 1: -m $opt{m}") if $opt{"m"} !~ /^\d+$/ or $opt{"m"} < 1;
	&exception::error("specify INT >= 0: -n $opt{n}") if $opt{"n"} !~ /^\d+$/;
	&exception::error("specify INT 1-" . main::max_depth . ": -a $opt{a}") if $opt{"a"} !~ /^\d+$/ or $opt{"a"} < 1 or $opt{"a"} > main::max_depth;
	&exception::error("specify FLOAT 0-1: -q $opt{q}") if !Scalar::Util::looks_like_number($opt{"q"}) or $opt{"q"} < 0 or $opt{"q"} > 1;
	$opt{"c"} = lc($opt{"c"}) if $opt{"c"};
	
	# 使用パラメータを表示
	print STDERR "Cutoff depth (cluster size): ", defined($opt{"d"}) ? $opt{"d"} : "recommended value","\n";
	print STDERR "Cutoff false discovery rate: $opt{q}\n";
	print STDERR "Use $opt{c} alignment (", $opt{"g"} > $opt{"e"} ? "affine" : "linear", " gap model)\n" if defined($opt{"c"});
	print STDERR "Gap opening penalty: $opt{g}\n" if defined($opt{"c"});
	print STDERR "Gap extension penalty: $opt{e}\n" if defined($opt{"c"});
	print STDERR "Match award: $opt{m}\n" if defined($opt{"c"});
	print STDERR "Mismatch penalty: $opt{n}\n" if defined($opt{"c"});
	print STDERR "Maximum amount of sequence reads to be aligned: $opt{a}\n" if defined($opt{"c"});
	
	# 入力ファイルのプレフィックスを取得
	my $input_prefix = shift(@ARGV);
	
	# 入力ファイルのプレフィックスを確認
	&exception::error("input file prefix not specified") unless defined($input_prefix);
	
	# 入力ファイルを確認
	&common::check_files(["$input_prefix.nsc", "$input_prefix.nso", "$input_prefix.nsr", "$input_prefix.nss"]);
	
	# アラインメント方法の対応番号を定義
	my %align_method = ("local" => 0, "global" => 1, "semi-global" => 2);
	
	# Nanoha Sequence Sketch (NSS) ファイルを開く
	open(NSS, "<", "$input_prefix.nss") or &exception::error("failed to open file: $input_prefix.nss");
	
	# NSSファイルをバイナリモードにする
	binmode(NSS);
	
	# NSSファイルのファイルヘッダーを読み込む
	read(NSS, my $nss_file_header, 8);
	
	# NSSファイルからパラメータ情報を読み込む
	read(NSS, my $parameter_info, 2);
	
	# ファイルバージョンの互換性を確認
	&exception::error("file version inconsistent with software version: $input_prefix.nss") if vec($parameter_info, 1, 4) >> 1 != [split(/\./, $version)]->[1] - 1;
	
	# ストランド特異性を取得
	my $strand_speicific = vec($parameter_info, 4, 1);
	
	# リード数を取得
	my $num_reads = ((-s "$input_prefix.nss") - 10) / (vec($parameter_info, 1, 8) + 2) / 4;
	
	# NSSファイルを閉じる
	close(NSS);
	
	# Nanoha Sequence Read (NSR) ファイルを開く
	open(NSR, "<", "$input_prefix.nsr") or &exception::error("failed to open file: $input_prefix.nsr");
	
	# NSRファイルをバイナリモードにする
	binmode(NSR);
	
	# NSRファイルのファイルヘッダーを読み込む
	read(NSR, my $nsr_file_header, 8);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsr, $input_prefix.nss") unless $nsr_file_header eq $nss_file_header;
	
	# シーケンスリードインデックスを取得
	seek(NSR, -($num_reads + 1) * 8, 2);
	read(NSR, my $seq_read_index, ($num_reads + 1) * 8);
	
	# NSRファイルを閉じる
	close(NSR);
	
	# Nanoha Sequence Orientation (NSO) ファイルを開く
	open(NSO, "<", "$input_prefix.nso") or &exception::error("failed to open file: $input_prefix.nso");
	
	# NSOファイルをバイナリモードにする
	binmode(NSO);
	
	# NSOファイルのファイルヘッダーを読み込む
	read(NSO, my $nso_file_header, 8);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nso, $input_prefix.nss") unless $nso_file_header eq $nss_file_header;
	
	# NSOファイルを読み込む
	print STDERR "Loading sequence orientations...";
	read(NSO, my $seq_ori, (-s "$input_prefix.nso") - 8);
	print STDERR "completed\n";
	
	# NSOファイルを閉じる
	close(NSO);
	
	# Nanoha Sequence Cluster (NSC) ファイルを開く
	open(NSC, "<", "$input_prefix.nsc") or &exception::error("failed to open file: $input_prefix.nsc");
	
	# NSCファイルをバイナリモードにする
	binmode(NSC);
	
	# NSCファイルのファイルヘッダーを読み込む
	read(NSC, my $nsc_file_header, 8);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsc, $input_prefix.nss") unless $nsc_file_header eq $nss_file_header;
	
	# NSCファイルからパラメータ情報を読み込む
	read(NSC, $parameter_info, 4);
	
	# 深度の閾値が未定義の場合はNSCファイルから読み込んだ値を使用
	$opt{"d"} = vec($parameter_info, 0, 32) unless defined($opt{"d"});
	
	# 変数を宣言
	my @worker_threads = ();
	my $num_error_threads = 0;
	
	# 入出力キューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 指定されたワーカースレッド数で並列処理
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_threads[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューを読み込みながら処理
			while (defined(my $serialized_data = $input->dequeue)) {
				# 変数を宣言
				my %read_seqs = ();
				
				# Nanoha Sequence Read (NSR) ファイルを開く
				open(NSR, "<", "$input_prefix.nsr") or &exception::error("failed to open file: $input_prefix.nsr");
				
				# NSRファイルをバイナリモードにする
				binmode(NSR);
				
				# 各リードについて処理
				foreach my $read_id (unpack("N*", $serialized_data)) {
					# NSRファイルのファイルポインタをシーケンスリードインデックスの位置にセット
					seek(NSR, vec($seq_read_index, $read_id - 1, 64), 0);
					
					# NSRファイルからリードシーケンスをバイナリ形式で読み込む
					read(NSR, my $read_seq, vec($seq_read_index, $read_id, 64) - tell(NSR));
					
					# リードシーケンスを復号
					$read_seq = &common::decode_sequence($read_seq);
					
					# 復号したシーケンスを追加 (-c指定時かつ逆鎖の場合は相補鎖に変換)
					$read_seqs{$read_id} = $opt{"c"} && vec($seq_ori, $read_id, 1) ? common::complementary($read_seq) : $read_seq;
				}
				
				# リード長の最小値と最大値を取得
				my ($shortest, $longest) = List::MoreUtils::minmax(map {length($_)} values(%read_seqs));
				
				# 順鎖リード数を算出
				my $num_plus_reads = grep {!vec($seq_ori, $_, 1)} keys(%read_seqs);
				
				# 逆鎖リード数を算出
				my $num_minus_reads = grep {vec($seq_ori, $_, 1)} keys(%read_seqs);
				
				# ストランドバイアス値を算出
				my $strand_bias = 1 - ($num_plus_reads < $num_minus_reads ? $num_plus_reads / $num_minus_reads : $num_minus_reads / $num_plus_reads);
				
				# リードシーケンスをリード長の降順にソート
				my @sorted_read_seqs = map {$read_seqs{$_}} sort {length($read_seqs{$b}) <=> length($read_seqs{$a}) || $a <=> $b} keys(%read_seqs);
				
				# リード長の長いものから指定した数のリードシーケンスを用いてコンセンサスリードシーケンスを構築 (-c指定時)
				@sorted_read_seqs = (&main::create_consensus_sequence(\@sorted_read_seqs, $opt{"a"}, $align_method{$opt{"c"}}, $opt{"m"}, -$opt{"n"}, -$opt{"g"}, -$opt{"e"})) if $opt{"c"};
				
				# リードシーケンスからギャップを除去
				s/-//g foreach @sorted_read_seqs;
				
				# 出力キューに最小リード長及び最大リード長、深度、ストランドバイアス値と各リードのリード長及びリードシーケンスをバイナリ形式で追加
				$output->enqueue(pack("LLLF(LA*)*", $shortest, $longest, scalar(keys(%read_seqs)), $strand_bias, map {length($_) => $_} @sorted_read_seqs));
			}
			
			# NSRファイルを閉じる
			close(NSR);
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 共有変数を宣言
	my $num_extracted_reads : shared;
	
	# データストリームスレッドを作成
	## ここからデータストリームスレッドの処理 ##
	$worker_threads[$opt{"p"}] = threads::async {
		# このスレッドのみを終了可能に変更
		threads->set_thread_exit_only(1);
		
		# 変数を宣言
		my $cluster_id = 0;
		
		# 出力キューを読み込みながら処理
		while (defined(my $serialized_data = $output->dequeue)) {
			# 変数を宣言
			my $read_id = 0;
			
			# クラスターIDを更新
			$cluster_id++;
			
			# 最小リード長及び最大リード長、深度、ストランドバイアス値と各リードのリードシーケンスを取得
			my ($shortest, $longest, $depth, $strand_bias, @read_seqs) = unpack("LLLF(L/A*)*", $serialized_data);
			
			# 各リードシーケンスをFASTA形式で出力
			++$read_id and print ">${input_prefix}_cluster_${cluster_id}", $opt{"c"} ? " shortest=$shortest longest=$longest depth=$depth strand_bias=$strand_bias" : ".$read_id", " length=", length($_), "\n$_\n" foreach @read_seqs;
			
			# 抽出したリード数を加算
			$num_extracted_reads += @read_seqs;
		}
		
		# スレッドを終了
		return(1);
	};
	## ここまでデータストリームスレッドの処理 ##
	
	# 変数を宣言
	my @clusters = ();
	
	# NSCファイルを読み込みながら処理
	print STDERR "Loading sequence clusters...";
	while (read(NSC, my $num_clustered_reads, 4)) {
		# NSCファイルから所属リード数×4バイト読み込んでクラスターリストに追加
		read(NSC, $clusters[scalar(@clusters)], vec($num_clustered_reads, 0, 32) * 4);
	}
	print STDERR "completed\n";
	
	# NSCファイルを閉じる
	close(NSC);
	
	# 変数を宣言
	my %binomial_coefficient = ();
	my @p_values = ();
	
	# 各クラスターについて処理
	print STDERR $opt{"c"} ? "Generating a consensus sequence read" : "Extracting sequence reads", " from each cluster...";
	foreach my $cluster (@clusters) {
		# ストランド特異性が有効な場合はp-valueに0を代入して以下の処理をスキップ
		push(@p_values, 0) and next if $strand_speicific;
		
		# 変数を宣言
		my @num_clustered_reads = (0, 0);
		
		# シーケンス方向ごとにリード数を集計
		$num_clustered_reads[vec($seq_ori, $_, 1)]++ foreach unpack("N*", $cluster);
		
		# リード数の合計を算出
		my $total_clustered_reads = List::Util::sum(@num_clustered_reads);
		
		# シーケンス方向で少ない方のリード数を取得
		my $less_clustered_reads = List::Util::min(@num_clustered_reads);
		
		# 二項係数を算出
		&calc_binomial_coefficient(\%binomial_coefficient, $total_clustered_reads, $less_clustered_reads);
		
		# ストランドバイアスのp-valueを算出
		push(@p_values, List::Util::sum(map {$binomial_coefficient{$total_clustered_reads}->[$_]} 0..$less_clustered_reads) / 2 ** $total_clustered_reads);
	}
	
	# 変数を宣言
	my @q_values = ();
	my @processed_clusters = ();
	
	# 各クラスターについてストランドバイアスのq-valueを算出
	push(@processed_clusters, $_) and $q_values[$_] = $p_values[$_] * @clusters / @processed_clusters foreach sort {$p_values[$a] <=> $p_values[$b] || $a <=> $b} 0..$#clusters;
	@q_values[reverse(@processed_clusters)] = map {$q_values[$_]} List::Util::reductions {$q_values[$a] < $q_values[$b] ? $a : $b} reverse(@processed_clusters);
	
	# ストランドバイアスのq-valueが指定値以上かつクラスターサイズが指定値以上の対象クラスターを取得
	my @target_clusters = grep {length($clusters[$_]) / 4 >= $opt{"d"}} grep {$strand_speicific or $q_values[$_] >= $opt{"q"}} 0..$#clusters;
	
	# 各対象クラスターについて実行中のスレッド数が指定値+1と一致している場合は所属リードのリードIDをバイナリ形式のまま入力キューに追加
	threads->list(threads::running) == $opt{"p"} + 1 and $input->enqueue($clusters[$_]) foreach @target_clusters;
	
	# 入力キューを終了
	$input->end;
	
	# 各ワーカースレッドが終了するまで待機
	$worker_threads[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 出力キューを終了
	$output->end;
	
	# データストリームスレッドが終了するまで待機
	$worker_threads[$opt{"p"}]->join or &exception::error("data stream thread abnormally exited");
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	# 抽出したリード数を表示
	print STDERR $num_extracted_reads, $opt{"c"} ? " consensus" : "", " sequence read", $num_extracted_reads > 1 ? "s" : "", $opt{"c"} ? " generated\n" : " extracted\n";
	
	return(1);
}

# 二項係数を算出 unify::calc_binomial_coefficient(二項係数ハッシュリファレンス, 整数, 整数)
sub calc_binomial_coefficient {
	# 引数を取得
	my ($binomial_coefficient, $n, $r) = @_;
	
	# 初期値を定義
	$binomial_coefficient->{$n}->[0] = Math::BigFloat->new(1) unless exists($binomial_coefficient->{$n});
	
	# 未算出の部分を算出
	$binomial_coefficient->{$n}->[$_] = $binomial_coefficient->{$n}->[$_ - 1] * ($n - $_ + 1) / $_ foreach scalar(@{$binomial_coefficient->{$n}})..$r;
	
	return(1);
}

## ここからconvertコマンドのパッケージ ##
package convert;

# コマンドとオプションを定義
sub define {
	$note = "Convert sequence reads from FASTA format to FASTQ format under specified conditions.";
	$usage = "<STDIN|in1.fa> [in2.fa ...] [>out.fq]";
	$option{"w"} = "\tUse 2-byte line feed code (CR+LF) for input files";
	return(1);
}

# コマンド本体
sub body {
	# 入力ファイルを確認
	&exception::error("input file not specified") unless @ARGV or -p STDIN;
	&common::check_files(\@ARGV);
	
	# 変数を宣言
	my $read_id = undef;
	my $read_seq = "";
	
	# 入力の改行コードを一時的に変更 (-w指定時)
	local $/ = "\r\n" if $opt{"w"};
	
	# FASTAファイルを読み込みながら処理
	print STDERR "Converting sequence reads from FASTA format to FASTQ format...";
	while (my $line = <>) {
		# 改行コードを除去
		chomp($line);
		
		# ヘッダー行でない場合はシーケンスを追加
		$read_seq .= $line if $line !~ /^>/;
		
		# リードIDが未定義でなくヘッダー行またはファイル末に到達した場合はリードシーケンスがnullでないか確認しFASTQ形式で出力してリードIDとともにリセット
		$read_seq and print join("\n", "@" . $read_id, $read_seq, "+", "I" x length($read_seq)), "\n" and ($read_id, $read_seq) = (undef, "") or &exception::error("null sequence found", $ARGV) if defined($read_id) and ($line =~ /^>/ or eof);
		
		# ヘッダー行の場合はリードIDを取得してnullでないか確認
		$read_id = substr($line, 1) or &exception::error("null header found", $ARGV) if $line =~ /^>/;
		
		# ファイル末に達した場合は読み込み行数をリセット
		$. = 0;
	}
	print STDERR "completed\n";
	
	return(1);
}

## ここからdumpコマンドのパッケージ ##
package dump;

# コマンドとオプションを定義
sub define {
	$note = "Dump k-mer minimizers counts under specified conditions.";
	$usage = "<prefix1> [prefix2 ...] [>out.tsv]";
	$option{"l STR "} = "Comma-separated sample label list";
	$option{"n INT "} = "Use n-th k-mer minimizers <1-> [1]";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("specify INT >= 1: -n $opt{n}") if ($opt{"n"} !~ /^\d+$/ or $opt{"n"} < 1);
	
	# 入力ファイルのプレフィックスを取得
	my @input_prefixes = @ARGV;
	
	# 入力ファイルのプレフィックスを確認
	&exception::error("input file prefix not specified") unless @input_prefixes;
	
	# 入力ファイルを確認
	&common::check_files([map {"$_.nss"} @input_prefixes]);
	
	# サンプルラベルを作成
	my @sample_labels = split(/,/, $opt{"l"}) if defined($opt{"l"});
	push(@sample_labels, $input_prefixes[$_]) foreach @sample_labels..@input_prefixes - 1;
	unshift(@sample_labels, "");
	
	# 変数を宣言
	my %minimizer_counts = ();
	my $kmer_size = undef;
	my $offset = ($opt{"n"} - 1) * 4;
	my $seq_read_index = "";
	
	# 各入力ファイルのプレフィックスについて処理
	for (my $i = 0;$i < @input_prefixes;$i++) {
		# Nanoha Sequence Sketch (NSS) ファイルを開く
		open(NSS, "<", "$input_prefixes[$i].nss") or &exception::error("failed to open file: $input_prefixes[$i].nss");
		
		# NSSファイルをバイナリモードにする
		binmode(NSS);
		
		# NSSファイルのファイルヘッダーを読み込む
		read(NSS, my $nss_file_header, 8);
		
		# NSSファイルからパラメータ情報を読み込む
		read(NSS, my $parameter_info, 2);
		
		# ファイルバージョンの互換性を確認
		&exception::error("file version inconsistent with software version: $input_prefixes[$i].nss") if vec($parameter_info, 1, 4) >> 1 != [split(/\./, $version)]->[1] - 1;
		
		# 未登録の場合はk-merサイズを登録
		$kmer_size = vec($parameter_info, 0, 4) unless $kmer_size;
		
		# k-merサイズが一致していることを確認
		&exception::error("k-mer size inconsistent: $input_prefixes[$i].nss") unless $kmer_size == vec($parameter_info, 0, 4);
		
		# minimizer数を取得
		my $num_minimizers = vec($parameter_info, 1, 8);
		
		# 指定値がminimizer数を超えていないことを確認
		&exception::error("number of minimizers less than $opt{n}: $input_prefixes[$i].nss") if $opt{"n"} > $num_minimizers;
		
		# リード数を取得
		my $num_reads = ((-s "$input_prefixes[$i].nss") - 10) / ($num_minimizers + 2) / 4;
		
		# シーケンススケッチインデックスを取得
		seek(NSS, -($num_reads + 1) * 8, 2);
		read(NSS, my $seq_sketch_index, ($num_reads + 1) * 8);
		
		# 各リードについて処理
		print STDERR "Loading sequence sketches...";
		for (my $read_id = 1;$read_id <= $num_reads;$read_id++) {
			# NSSファイルのファイルポインタをシーケンススケッチインデックス及びオフセット値の示す位置にセット
			seek(NSS, vec($seq_sketch_index, $read_id, 64) + $offset, 0);
			
			# NSSファイルからminimizerを読み込む
			read(NSS, my $minimizer, 4);
			
			# minimizerが未登録の場合は登録
			$minimizer_counts{$minimizer} = pack("N*", (0) x @input_prefixes) unless exists($minimizer_counts{$minimizer});
			
			# minimizer出現回数のカウントを加算
			vec($minimizer_counts{$minimizer}, $i, 32)++;
		}
		print STDERR "completed\n";
		
		# NSSファイルを閉じる
		close(NSS);
	}
	
	# サンプルラベルを出力
	print join("\t", @sample_labels), "\n";
	
	# 各minimizerについて処理
	foreach my $minimizer (sort(keys(%minimizer_counts))) {
		# minimizerの配列を取得
		my $minimizer_seq = substr($convert->to_base(vec($minimizer, 0, 32)), 1);
		
		# 各サンプルのminimizerカウントを取得
		my @minimizer_count = unpack("N*", $minimizer_counts{$minimizer});
		
		# 各サンプルのminimizerカウントを出力
		print join("\t", $minimizer_seq, @minimizer_count), "\n";
	}
	
	return(1);
}

## ここから共通処理のパッケージ ##
package common;

# ファイル確認 common::check_files(ファイルリストリファレンス)
sub check_files {
	# 引数を取得
	my ($files) = @_;
	
	# 各ファイルについて処理
	foreach my $file (grep {defined($_)} @{$files}) {
		&exception::error("file not found: $file") unless -f $file;
		&exception::error("file unreadable: $file") unless -r $file;
		&exception::error("null file specified: $file") unless -s $file;
	}
	return(1);
}

# 相補鎖変換 common::complementary(配列)
sub complementary {
	# 引数を取得
	my ($seq) = @_;
	
	# シーケンスを逆順に並べ替える
	my $complementary_seq = reverse($seq);
	
	# 相補的な塩基に置換
	$complementary_seq =~ tr/ATGCRYKMDBVH/TACGYRMKHVBD/;
	
	# 相補鎖シーケンスを返す
	return($complementary_seq);
}

# シーケンスを復号 common::decode_sequence(配列)
sub decode_sequence {
	# 引数を取得
	my ($seq) = @_;
	
	# 復号したシーケンスを返す
	return(join("", map {$convert->to_base(vec($seq, $_ + 16, 2))} 0..vec($seq, 0, 32) - 1));
}

# サブルーチンを追加

# パッケージを追加
### 編集範囲 終了 ###
__END__
__C__
// マクロを定義
#define newAV_mortal sv_2mortal(newAV())

// LLCSに基づくJaccard類似度を算出 main::calc_LLCS_Jaccard_similarities(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性)
AV* calc_LLCS_Jaccard_similarities(AV* read_seqs, AV* read_lens, unsigned int strand_speicific) {
	// 戻り値のリストを作成
	AV* LLCS_Jaccard_similarities = newAV_mortal;
	
	// クエリーリードシーケンスを取得
	const uint64_t* query_seq = SvPV_nolen(*av_fetch(read_seqs, 0, 0));
	
	// クエリーリード長を取得
	const uint32_t query_len = SvUV(*av_fetch(read_lens, 0, 0));
	
	// クエリーブロック数を算出
	const uint32_t num_query_blocks = (query_len >> 6) + ((query_len & 0x3F) > 0);
	
	// 変数を宣言
	uint64_t* query_matrix[4];
	uint64_t query_base[2];
	uint64_t llcs[2];
	uint64_t mask = 0xFFFFFFFFFFFFFFFF;
	uint64_t* v;
	uint64_t u;
	uint64_t w;
	
	// メモリを確保
	query_matrix[0] = malloc(sizeof(uint64_t) * num_query_blocks);
	query_matrix[1] = malloc(sizeof(uint64_t) * num_query_blocks);
	query_matrix[2] = malloc(sizeof(uint64_t) * num_query_blocks);
	query_matrix[3] = malloc(sizeof(uint64_t) * num_query_blocks);
	v = malloc(sizeof(uint64_t) * num_query_blocks);
	
	// メモリが確保されたことを確認
	if (query_matrix[0] == NULL || query_matrix[1] == NULL || query_matrix[2] == NULL || query_matrix[3] == NULL || v == NULL) {goto CLEANUP;}
	
	// クエリー行列を作成
	for (uint32_t i = 0;i < 4;i++) {
		for (uint32_t j = 0;j < num_query_blocks;j++) {
			query_base[0] = query_seq[j * 2] ^ mask;
			query_base[0] &= query_base[0] >> 1;
			query_base[1] = query_seq[j * 2 + 1] ^ mask;
			query_base[1] &= query_base[1] << 1;
			query_matrix[i][j] = 0;
			for (uint32_t k = 0;k < 32;k++) {
				query_matrix[i][j] |= query_base[0] >> k & 0x0000000000000001 << k | query_base[1] << k & 0x8000000000000000 >> k;
			}
		}
		mask -= 0x5555555555555555;
	}
	
	// クエリー行列の端をマスク
	query_matrix[0][num_query_blocks - 1] &= 0xFFFFFFFFFFFFFFFF >> (-query_len & 0x3F);
	
	// 残りの各リードについて処理
	for (uint32_t n = 1;n <= av_len(read_seqs);n++) {
		// リードシーケンスを取得
		uint64_t* read_seq = SvPV_nolen(*av_fetch(read_seqs, n, 0));
		
		// リード長を取得
		uint32_t read_len = SvUV(*av_fetch(read_lens, n, 0));
		
		// ブロック数を算出
		uint32_t num_blocks = (read_len >> 5) + ((read_len & 0x1F) > 0);
		
		// 順鎖及び逆鎖について処理
		for (uint32_t strand = 0;strand < 2;strand++) {
			// 計算要素を初期化
			llcs[strand] = 0;
			for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = 0xFFFFFFFFFFFFFFFF;}
			
			// LLCSを算出
			for (uint32_t i = strand * (num_blocks - 1);i < num_blocks;i += 1 - 2 * strand) {
				uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
				for (uint32_t j = strand * (end - 2);j < end;j += 2 - 4 * strand) {
					uint32_t base = read_seq[i] >> j & 0x03;
					uint64_t carry = 0;
					for (uint32_t k = 0;k < num_query_blocks;k++) {
						u = v[k] & query_matrix[base][k];
						w = v[k] + u + carry;
						carry = w < v[k];
						v[k] = w | v[k] - u;
					}
				}
			}
			for (uint32_t i = 0;i < num_query_blocks;i++) {
				w = ~v[i];
				w = (w & 0x5555555555555555) + (w >> 1 & 0x5555555555555555);
				w = (w & 0x3333333333333333) + (w >> 2 & 0x3333333333333333);
				w = (w & 0x0F0F0F0F0F0F0F0F) + (w >> 4 & 0x0F0F0F0F0F0F0F0F);
				w = (w & 0x00FF00FF00FF00FF) + (w >> 8 & 0x00FF00FF00FF00FF);
				w = (w & 0x0000FFFF0000FFFF) + (w >> 16 & 0x0000FFFF0000FFFF);
				w = (w & 0x00000000FFFFFFFF) + (w >> 32 & 0x00000000FFFFFFFF);
				llcs[strand] += w;
			}
			
			// ストランド特異性が有効な場合はループを終了
			if (strand_speicific) {break;}
			
			// リードシーケンスを相補的な塩基に変換
			for (uint32_t i = 0;i < num_blocks;i++) {read_seq[i] = ~read_seq[i];}
			
		}
		// ストランド特異性が有効な場合は順鎖に対するLLCSに基づくJaccard類似度を算出して戻り値に追加し以下の処理をスキップ
		if (strand_speicific) {
			av_push(LLCS_Jaccard_similarities, newSVnv((double)llcs[0] / (query_len + read_len - llcs[0])));
			continue;
		}
		
		// LLCSに基づくJaccard類似度を算出して戻り値に追加 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		av_push(LLCS_Jaccard_similarities, newSVnv((double)(1 - 2 * (llcs[0] < llcs[1])) * llcs[llcs[0] < llcs[1]] / (query_len + read_len - llcs[llcs[0] < llcs[1]])));
	}
	
	CLEANUP:
	// メモリを解放
	free(query_matrix[0]);
	free(query_matrix[1]);
	free(query_matrix[2]);
	free(query_matrix[3]);
	free(v);
	
	// LLCSに基づくJaccard類似度を返す
	return LLCS_Jaccard_similarities;
}

__C__
// 共通のマクロを定義
#define newAV_mortal sv_2mortal(newAV())

// SSE4.1用のマクロ及び大域変数を定義
#if word_size_index == 1
	#define vector_int __m128i
	#define vector_double __m128d
	#define vec_zero _mm_setzero_si128()
	#define vec_load(a) _mm_load_si128(a)
	#define vec_cast_dword(a) _mm_cvtsi128_si32(a)
	#define vec_set_qword(a, b) _mm_set_epi64x(a, b)
	#define vec_splat_qword(a) _mm_set1_epi64x(a)
	#define vec_splat_byte(a) _mm_set1_epi8(a)
	#define vec_movemask_double(a) _mm_movemask_pd(a)
	#define vec_shuffle_dword(a, b) _mm_shuffle_epi32(a, b)
	#define vec_shuffle_byte(a, b) _mm_shuffle_epi8(a, b)
	#define vec_pack_word(a, b) _mm_packus_epi16(a, b)
	#define vec_cmpeq_qword(a, b) _mm_cmpeq_epi64(a, b)
	#define vec_cmpgt_dword(a, b) _mm_cmpgt_epi32(a, b)
	#define vec_add_qword(a, b) _mm_add_epi64(a, b)
	#define vec_add_dword(a, b) _mm_add_epi32(a, b)
	#define vec_add_byte(a, b) _mm_add_epi8(a, b)
	#define vec_sub_qword(a, b) _mm_sub_epi64(a, b)
	#define vec_sub_dword(a, b) _mm_sub_epi32(a, b)
	#define vec_mad_word(a, b) _mm_maddubs_epi16(a, b)
	#define vec_avg_byte(a, b) _mm_avg_epu8(a, b)
	#define vec_sad(a, b) _mm_sad_epu8(a, b)
	#define vec_srli_qword(a, b) _mm_srli_epi64(a, b)
	#define vec_srai_dword(a, b) _mm_srai_epi32(a, b)
	#define vec_and(a, b) _mm_and_si128(a, b)
	#define vec_andnot(a, b) _mm_andnot_si128(a, b)
	#define vec_or(a, b) _mm_or_si128(a, b)
	#define vec_xor(a, b) _mm_xor_si128(a, b)
	
// AVX2用のマクロ及び大域変数を定義
#elif word_size_index == 2
	#define vector_int __m256i
	#define vector_double __m256d
	#define vec_zero _mm256_setzero_si256()
	#define vec_load(a) _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_load_si128(a)), _mm_load_si128(a + 4), 1)
	#define vec_cast_dword(a) _mm256_cvtsi256_si32(a)
	#define vec_set_qword(a, b, c, d) _mm256_set_epi64x(a, b, c, d)
	#define vec_splat_qword(a) _mm256_set1_epi64x(a)
	#define vec_splat_byte(a) _mm256_set1_epi8(a)
	#define vec_movemask_double(a) _mm256_movemask_pd(a)
	#define vec_permute_dword(a, b) _mm256_permutevar8x32_epi32(a, b)
	#define vec_shuffle_dword(a, b) _mm256_shuffle_epi32(a, b)
	#define vec_shuffle_byte(a, b) _mm256_shuffle_epi8(a, b)
	#define vec_pack_word(a, b) _mm256_packus_epi16(a, b)
	#define vec_cmpeq_qword(a, b) _mm256_cmpeq_epi64(a, b)
	#define vec_cmpgt_dword(a, b) _mm256_cmpgt_epi32(a, b)
	#define vec_add_qword(a, b) _mm256_add_epi64(a, b)
	#define vec_add_dword(a, b) _mm256_add_epi32(a, b)
	#define vec_add_byte(a, b) _mm256_add_epi8(a, b)
	#define vec_sub_qword(a, b) _mm256_sub_epi64(a, b)
	#define vec_sub_dword(a, b) _mm256_sub_epi32(a, b)
	#define vec_mad_word(a, b) _mm256_maddubs_epi16(a, b)
	#define vec_avg_byte(a, b) _mm256_avg_epu8(a, b)
	#define vec_sad(a, b) _mm256_sad_epu8(a, b)
	#define vec_srli_qword(a, b) _mm256_srli_epi64(a, b)
	#define vec_srai_dword(a, b) _mm256_srai_epi32(a, b)
	#define vec_and(a, b) _mm256_and_si256(a, b)
	#define vec_andnot(a, b) _mm256_andnot_si256(a, b)
	#define vec_or(a, b) _mm256_or_si256(a, b)
	#define vec_xor(a, b) _mm256_xor_si256(a, b)
	
// AVX-512用のマクロ及び大域変数を定義
#elif word_size_index == 3
	#define vector_int __m512i
	#define vector_double __m512d
	#define vec_zero _mm512_setzero_si512()
	#define vec_load(a) _mm512_inserti32x4(_mm512_inserti32x4(_mm512_inserti32x4(_mm512_castsi128_si512(_mm_load_si128(a)), _mm_load_si128(a + 4), 1), _mm_load_si128(a + 8), 2), _mm_load_si128(a + 12), 3)
	#define vec_cast_dword(a) _mm512_cvtsi512_si32(a)
	#define vec_set_qword(a, b, c, d, e, f, g, h) _mm512_set_epi64(a, b, c, d, e, f, g, h)
	#define vec_splat_qword(a) _mm512_set1_epi64(a)
	#define vec_splat_byte(a) _mm512_set1_epi8(a)
	#define vec_permute_qword(a, b) _mm512_permutex_epi64(a, b)
	#define vec_permute_dword(a, b) _mm512_permutexvar_epi32(b, a)
	#define vec_shuffle_dword(a, b) _mm512_shuffle_epi32(a, b)
	#define vec_shuffle_byte(a, b) _mm512_shuffle_epi8(a, b)
	#define vec_pack_word(a, b) _mm512_packus_epi16(a, b)
	#define vec_cmpeq_qword(a, b) _mm512_cmpeq_epu64_mask(a, b)
	#define vec_cmplt_qword(a, b) _mm512_cmplt_epu64_mask(a, b)
	#define vec_cmpgt_dword(a, b) _mm512_cmpgt_epi32_mask(a, b)
	#define vec_add_qword(a, b) _mm512_add_epi64(a, b)
	#define vec_add_dword(a, b) _mm512_add_epi32(a, b)
	#define vec_add_byte(a, b) _mm512_add_epi8(a, b)
	#define vec_sub_qword(a, b) _mm512_sub_epi64(a, b)
	#define vec_sub_dword(a, b) _mm512_sub_epi32(a, b)
	#define vec_mask_sub_qword(a, b, c, d) _mm512_mask_sub_epi64(a, b, c, d)
	#define vec_mask_sub_dword(a, b, c, d) _mm512_mask_sub_epi32(a, b, c, d)
	#define vec_mad_word(a, b) _mm512_maddubs_epi16(a, b)
	#define vec_avg_byte(a, b) _mm512_avg_epu8(a, b)
	#define vec_sad(a, b) _mm512_sad_epu8(a, b)
	#define vec_srli_qword(a, b) _mm512_srli_epi64(a, b)
	#define vec_srai_dword(a, b) _mm512_srai_epi32(a, b)
	#define vec_and(a, b) _mm512_and_si512(a, b)
	#define vec_andnot(a, b) _mm512_andnot_si512(a, b)
	#define vec_or(a, b) _mm512_or_si512(a, b)
	#define vec_xor(a, b) _mm512_xor_si512(a, b)
	#define vec_ternarylogic(a, b, c, d) _mm512_ternarylogic_epi64(a, b, c, d)
#endif

// LLCSに基づくJaccard類似度をSIMD演算を用いて算出 main::calc_LLCS_Jaccard_similarities_by_simd(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性)
AV* calc_LLCS_Jaccard_similarities_by_simd(AV* read_seqs, AV* read_lens, unsigned int strand_speicific) {
	// 戻り値のリストを作成
	AV* LLCS_Jaccard_similarities = newAV_mortal;
	
	// クエリーリードシーケンスを取得
	const uint64_t* query_seq = SvPV_nolen(*av_fetch(read_seqs, 0, 0));
	
	// クエリーリード長を取得
	const uint32_t query_len = SvUV(*av_fetch(read_lens, 0, 0));
	
	// クエリーブロック数を算出
	const uint32_t num_query_blocks = (query_len >> word_size_index + 6) + ((query_len & ((0x40 << word_size_index) - 1)) > 0);
	
	// 定数を定義
	#if word_size_index == 1
		const vector_int compress_table1 = vec_set_qword(0x0F0E0B0A0D0C0908, 0x0706030205040100);
		const vector_int compress_table2 = vec_set_qword(0x1001100110011001, 0x1001100110011001);
		const vector_int maskend_table1 = vec_set_qword(0xF1F2F3F4F5F6F7F8, 0xF9FAFBFCFDFEFF00);
		const vector_int maskend_table2 = vec_set_qword(0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF);
		const vector_int extract_table = vec_set_qword(0x0000000000000002, 0x0000000000000001);
		const vector_int popcnt_table1 = vec_set_qword(0x0403030203020201, 0x0302020102010100);
		const vector_int popcnt_table2 = vec_set_qword(0x0405050605060607, 0x0506060706070708);
	#elif word_size_index == 2
		const vector_int compress_table1 = vec_set_qword(0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100);
		const vector_int compress_table2 = vec_set_qword(0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001);
		const vector_int maskend_table1 = vec_set_qword(0xE1E2E3E4E5E6E7E8, 0xE9EAEBECEDEEEFF0, 0xF1F2F3F4F5F6F7F8, 0xF9FAFBFCFDFEFF00);
		const vector_int maskend_table2 = vec_set_qword(0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF);
		const vector_int extract_table = vec_set_qword(0x0000000000000008, 0x0000000000000004, 0x0000000000000002, 0x0000000000000001);
		const vector_int popcnt_table1 = vec_set_qword(0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100);
		const vector_int popcnt_table2 = vec_set_qword(0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708);
		const vector_int permute_table = vec_set_qword(0x0000000300000002, 0x0000000100000000, 0x0000000700000006, 0x0000000500000004);
	#elif word_size_index == 3
		const vector_int compress_table1 = vec_set_qword(0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100);
		const vector_int compress_table2 = vec_set_qword(0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001);
		const vector_int maskend_table1 = vec_set_qword(0xC1C2C3C4C5C6C7C8, 0xC9CACBCCCDCECFD0, 0xD1D2D3D4D5D6D7D8, 0xD9DADBDCDDDEDFE0, 0xE1E2E3E4E5E6E7E8, 0xE9EAEBECEDEEEFF0, 0xF1F2F3F4F5F6F7F8, 0xF9FAFBFCFDFEFF00);
		const vector_int maskend_table2 = vec_set_qword(0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF);
		const vector_int popcnt_table1 = vec_set_qword(0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100);
		const vector_int popcnt_table2 = vec_set_qword(0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708);
		const vector_int permute_table = vec_set_qword(0x0000000700000006, 0x0000000500000004, 0x0000000300000002, 0x0000000100000000, 0x0000000F0000000E, 0x0000000D0000000C, 0x0000000B0000000A, 0x0000000900000008);
	#endif
	const vector_int mask_ff = vec_splat_qword(0xFFFFFFFFFFFFFFFF);
	const vector_int mask_55 = vec_splat_qword(0x5555555555555555);
	const vector_int mask_0f = vec_splat_qword(0x0F0F0F0F0F0F0F0F);
	
	// 変数を宣言
	vector_int* query_matrix[4];
	vector_int query_base[2];
	vector_int llcs[2];
	vector_int mask = mask_ff;
	vector_int* v;
	vector_int u;
	vector_int w;
	vector_int x;
	int32_t signed_llcs;
	
	// メモリを確保
	query_matrix[0] = _mm_malloc(sizeof(vector_int) * num_query_blocks, sizeof(vector_int));
	query_matrix[1] = _mm_malloc(sizeof(vector_int) * num_query_blocks, sizeof(vector_int));
	query_matrix[2] = _mm_malloc(sizeof(vector_int) * num_query_blocks, sizeof(vector_int));
	query_matrix[3] = _mm_malloc(sizeof(vector_int) * num_query_blocks, sizeof(vector_int));
	v = _mm_malloc(sizeof(vector_int) * num_query_blocks, sizeof(vector_int));
	
	// メモリが確保されたことを確認
	if (query_matrix[0] == NULL || query_matrix[1] == NULL || query_matrix[2] == NULL || query_matrix[3] == NULL || v == NULL) {goto CLEANUP;}
	
	// クエリー行列を作成
	for (uint32_t i = 0;i < 4;i++) {
		for (uint32_t j = 0;j < num_query_blocks;j++) {
			for (uint32_t k = 0;k < 2;k++) {
				query_base[k] = vec_load(query_seq + (j * sizeof(vector_int) * 2 + k * sizeof(__m128i)) / sizeof(query_seq));
				query_base[k] = vec_xor(query_base[k], mask);
				query_base[k] = vec_and(vec_and(vec_srli_qword(query_base[k], 1), query_base[k]), mask_55);
				query_base[k] = vec_and(vec_or(vec_srli_qword(query_base[k], 3), query_base[k]), mask_0f);
				query_base[k] = vec_mad_word(vec_shuffle_byte(compress_table1, query_base[k]), compress_table2);
			}
			query_matrix[i][j] = vec_pack_word(query_base[0], query_base[1]);
		}
		mask = vec_sub_qword(mask, mask_55);
	}
	
	// クエリー行列の端をマスク
	mask = vec_add_byte(maskend_table1, vec_splat_byte(query_len - 1 >> 3 & ((0x08 << word_size_index) - 1)));
	mask = vec_or(vec_and(vec_srli_qword(vec_avg_byte(vec_xor(mask, mask_ff), vec_zero), 4), mask_0f), vec_andnot(mask_0f, mask));
	mask = vec_shuffle_byte(maskend_table2, vec_or(mask, vec_splat_byte(query_len & 0x07)));
	query_matrix[0][num_query_blocks - 1] = vec_and(query_matrix[0][num_query_blocks - 1], mask);
	
	// 残りの各リードについて処理
	for (uint32_t n = 1;n <= av_len(read_seqs);n++) {
		// リードシーケンスを取得
		uint64_t* read_seq = SvPV_nolen(*av_fetch(read_seqs, n, 0));
		
		// リード長を取得
		uint32_t read_len = SvUV(*av_fetch(read_lens, n, 0));
		
		// ブロック数を算出
		uint32_t num_blocks = (read_len >> 5) + ((read_len & 0x1F) > 0);
		
		// 順鎖及び逆鎖について処理
		for (uint32_t strand = 0;strand < 2;strand++) {
			// 計算要素を初期化
			llcs[strand] = vec_zero;
			for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = mask_ff;}
			
			// LLCSを算出
			for (uint32_t i = strand * (num_blocks - 1);i < num_blocks;i += 1 - 2 * strand) {
				uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
				for (uint32_t j = strand * (end - 2);j < end;j += 2 - 4 * strand) {
					uint32_t base = read_seq[i] >> j & 0x03;
					uint64_t carry = 0;
					for (uint32_t k = 0;k < num_query_blocks;k++) {
						u = vec_and(v[k], query_matrix[base][k]);
						w = vec_add_qword(v[k], u);
						#if word_size_index <= 2
							carry = vec_movemask_double((vector_double)vec_cmpeq_qword(w, mask_ff)) + (vec_movemask_double((vector_double)vec_or(vec_andnot(w, v[k]), u)) << 1) + (carry >> (0x01 << word_size_index));
							v[k] = vec_or(vec_sub_qword(w, vec_cmpeq_qword(vec_and(vec_splat_qword(carry), extract_table), extract_table)), vec_xor(v[k], u));
						#else
							carry = vec_cmpeq_qword(w, mask_ff) + (vec_cmplt_qword(w, v[k]) << 1) + (carry >> (0x01 << word_size_index));
							v[k] = vec_or(vec_mask_sub_qword(w, carry, w, mask_ff), vec_xor(v[k], u));
						#endif
					}
				}
			}
			u = num_query_blocks & 0x01 ? v[0] : mask_ff;
			for (uint32_t i = num_query_blocks & 0x01;i < num_query_blocks;i += 2) {
				#if word_size_index <= 2
					x = vec_xor(v[i], v[i + 1]);
					w = vec_or(vec_and(v[i], v[i + 1]), vec_and(u, x));
					u = vec_xor(u, x);
				#else
					w = vec_ternarylogic(u, v[i], v[i + 1], 0xE8);
					u = vec_ternarylogic(u, v[i], v[i + 1], 0x96);
				#endif
				llcs[strand] = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(w, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(w, 4), mask_0f))), llcs[strand]);
			}
			llcs[strand] = vec_add_qword(llcs[strand], llcs[strand]);
			llcs[strand] = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(u, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(u, 4), mask_0f))), llcs[strand]);
			
			// ストランド特異性が有効な場合はループを終了
			if (strand_speicific) {break;}
			
			// リードシーケンスを相補的な塩基に変換
			for (uint32_t i = 0;i < num_blocks;i++) {read_seq[i] = ~read_seq[i];}
		}
		
		// ストランド特異性が有効な場合は順鎖に対するLLCSに基づくJaccard類似度を算出して戻り値に追加し以下の処理をスキップ
		if (strand_speicific) {
			llcs[0] = vec_add_dword(vec_shuffle_dword(llcs[0], 0x4E), llcs[0]);
			#if word_size_index >= 3
				llcs[0] = vec_add_dword(vec_permute_qword(llcs[0], 0x4E), llcs[0]);
			#endif
			#if word_size_index >= 2
				llcs[0] = vec_add_dword(vec_permute_dword(llcs[0], permute_table), llcs[0]);
			#endif
			signed_llcs = vec_cast_dword(llcs[0]);
			av_push(LLCS_Jaccard_similarities, newSVnv((double)signed_llcs / (query_len + read_len - signed_llcs)));
			continue;
		}
		
		// LLCSに基づくJaccard類似度を算出して戻り値に追加 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		w = vec_add_dword(llcs[0], llcs[1]);
		w = vec_add_dword(vec_shuffle_dword(w, 0x4E), w);
		#if word_size_index >= 3
			w = vec_add_dword(vec_permute_qword(w, 0x4E), w);
		#endif
		#if word_size_index >= 2
			w = vec_add_dword(vec_permute_dword(w, permute_table), w);
		#endif
		u = vec_sub_dword(llcs[0], llcs[1]);
		u = vec_add_dword(vec_shuffle_dword(u, 0x4E), u);
		#if word_size_index >= 3
			u = vec_add_dword(vec_permute_qword(u, 0x4E), u);
		#endif
		#if word_size_index >= 2
			u = vec_add_dword(vec_permute_dword(u, permute_table), u);
		#endif
		#if word_size_index <= 2
			x = vec_cmpgt_dword(vec_zero, u);
			w = vec_sub_dword(vec_xor(w, x), x);
		#else
			w = vec_mask_sub_dword(w, vec_cmpgt_dword(vec_zero, u), vec_xor(w, mask_ff), mask_ff);
		#endif
		signed_llcs = vec_cast_dword(vec_srai_dword(vec_add_dword(w, u), 1));
		av_push(LLCS_Jaccard_similarities, newSVnv((double)signed_llcs / (query_len + read_len - abs(signed_llcs))));
	}
	
	CLEANUP:
	// メモリを解放
	_mm_free(query_matrix[0]);
	_mm_free(query_matrix[1]);
	_mm_free(query_matrix[2]);
	_mm_free(query_matrix[3]);
	_mm_free(v);
	
	// LLCSに基づくJaccard類似度を返す
	return LLCS_Jaccard_similarities;
}

// コンパイル時に設定したワードサイズインデックスを返す main::get_word_size_index()
unsigned int get_word_size_index() {return word_size_index;}

__CPP__
// マクロを定義
#define newAV_mortal (AV*)sv_2mortal((SV*)newAV())

// weighted minimizerを生成 main::generate_weighted_minimizers(基数変換したk-merリストリファレンス, minimizer生成個数)
AV* generate_weighted_minimizers(AV* converted_kmers, unsigned char num_minimizers) {
	// 戻り値のリストを作成
	AV* minimizers = newAV_mortal;
	
	// 定数を定義
	const uint32_t num_kmers = av_len(converted_kmers) + 1;
	
	// 変数を宣言
	std::unordered_map<uint32_t, std::vector<double>> hash_values;
	
	// 0以上1未満の実数値の一様分布器を作成
	std::uniform_real_distribution<> rnd(0.0, 1.0);
	
	// 基数変換した各k-merについて処理
	for (uint32_t i = 0;i < num_kmers;i += 2) {
		// 基数変換したk-merを取得
		uint32_t converted_kmer = SvUV(*av_fetch(converted_kmers, i, 0));
		
		// k-merの出現回数を取得
		double count = SvNV(*av_fetch(converted_kmers, i + 1, 0));
		
		// 指定した個数のハッシュ値を登録する配列を宣言
		hash_values[converted_kmer].resize(num_minimizers);
		
		// 基数変換したk-merをシード値に用いてメルセンヌ・ツイスターによる擬似乱数生成器を作成
		std::mt19937 mt(converted_kmer);
		
		// 指定した個数のハッシュ値を生成
		for (uint32_t j = 0;j < num_minimizers;j++) {hash_values[converted_kmer][j] = 1 - std::pow(rnd(mt), 1 / count);}
	}
	
	// 指定した個数分だけ処理
	for (uint32_t i = 0;i < num_minimizers;i++) {
		// minimizerを生成
		uint32_t minimizer = std::min_element(hash_values.begin(), hash_values.end(), [i](const auto& x, const auto& y) {return x.second[i] < y.second[i];})->first;
		
		// 生成したminimizerを戻り値に追加
		av_push(minimizers, newSVuv(minimizer));
	}
	
	// minimizerを返す
	return minimizers;
}

__CPP__
// マクロを定義
#define newAV_mortal (AV*)sv_2mortal((SV*)newAV())

// minimizer一致数ごとにminimizersをバケットソート main::sort_minimizers(minimizerバケットリストリファレンス, minimizers, minimizer一致数閾値)
AV* sort_minimizers(AV* minimizer_buckets, SV* minimizers_sv, unsigned char cutoff_num_matched_minimizers) {
	// 戻り値のリストを作成
	AV* sorted_minimizers = newAV_mortal;
	
	// minimizer数を取得
	uint8_t num_minimizers = av_len(minimizer_buckets) + 1;
	
	// minimizerを取得
	char* minimizers = SvPV_nolen(minimizers_sv);
	
	// 変数を宣言
	std::unordered_map<std::string, uint8_t> num_matched_minimizers;
	
	// minimizer一致数ごとにminimizersをバケットソート
	for (uint8_t i = 0;i < num_minimizers;i++) {
		// minimizerバケットを取得
		HV* minimizer_bucket = (HV*)SvRV(*av_fetch(minimizer_buckets, i, 0));
		
		// minimizerが一致するminimizersリストが存在しない場合は以下の処理をスキップ
		if (!hv_exists(minimizer_bucket, minimizers + i * 4, 4)) {continue;}
		
		// minimizerが一致するminimizersリストを取得
		SV* matched_minimizers_list_sv = *hv_fetch(minimizer_bucket, minimizers + i * 4, 4, 0);
		char* matched_minimizers_list = SvPV_nolen(matched_minimizers_list_sv);
		
		// minimizerが一致する各minimizersについて処理
		for (uint32_t j = 0;j < SvCUR(matched_minimizers_list_sv);j += num_minimizers * 4) {
			// minimizerが一致するminimizersを取得
			std::string matched_minimizers(matched_minimizers_list + j, num_minimizers * 4);
			
			// minimizer一致数を加算
			num_matched_minimizers[matched_minimizers]++;
		}
	}
	
	// 戻り値リストに空リストリファレンスを作成して追加
	for (uint8_t i = cutoff_num_matched_minimizers;i <= num_minimizers;i++) {av_push(sorted_minimizers, (SV*)newRV((SV*)newAV_mortal));}
	
	// 各minimizersのうち、指定した個数以上のminimizerが一致するものについてminimizer一致数に対応するバケットにminimizersを追加
	for (const auto& num_matched_minimizer : num_matched_minimizers) {
		if (num_matched_minimizer.second >= cutoff_num_matched_minimizers) {
			av_push((AV*)SvRV(*av_fetch(sorted_minimizers, num_matched_minimizer.second - cutoff_num_matched_minimizers, 0)), newSVpv(num_matched_minimizer.first.data(), num_matched_minimizer.first.size()));
		}
	}
	
	// 結果を返す
	return sorted_minimizers;
}

__CPP__
// 密度を算出 main::calc_density(内部エッジの重みの総和, ノード数)
double calc_density(double sum_internal_edges, unsigned int num_nodes) {
	// 密度を算出
	double density = num_nodes > 1 ? sum_internal_edges / num_nodes / (num_nodes - 1) : 0.0;
	
	// 密度を返す
	return density;
}

// 部分的なQxを算出 main::calc_partial_Qx(内部エッジの重みの総和, ノード数, クラスター次数, グラフ次数, グラフ密度)
double calc_partial_Qx(double sum_internal_edges, unsigned int num_nodes, double cluster_degree, double graph_degree, double graph_density) {
	// 密度指数を算出
	double density_index = calc_density(sum_internal_edges, num_nodes) - graph_density;
	
	// 部分的なQxを算出
	double partial_Qx = sum_internal_edges * density_index - pow(cluster_degree * density_index, 2) / graph_degree;
	
	// 部分的なQxを返す
	return partial_Qx;
}

// excess modularity density (Qx) に基づくリードのクラスタリング main::create_sequence_clusters(クラスター割り当てリストリファレンス, 極大クラスターリストリファレンス, シーケンスグラフインデックス, nanohaシーケンスグラフファイルパス)
double create_sequence_clusters(SV* cluster_assignment_sv, AV* max_clusters_av, SV* seq_graph_index_sv, char* nanoha_sequence_graph) {
	// クラスター割り当てリストを取得
	uint32_t* cluster_assignment = (uint32_t*)SvPV_nolen(cluster_assignment_sv);
	
	// リード数を取得
	uint32_t num_reads = SvCUR(cluster_assignment_sv) / sizeof(uint32_t) - 1;
	
	// シーケンスグラフインデックスを取得
	uint64_t* seq_graph_index = (uint64_t*)SvPV_nolen(seq_graph_index_sv);
	
	// NSGファイルを開く
	FILE* NSG = fopen(nanoha_sequence_graph, "rb");
	
	// ファイルが開けなかった場合はNaNを返す
	if (NSG == NULL) {return log(0.0) / log(0.0);}
	
	// 変数を宣言
	std::vector<std::vector<uint32_t>> max_clusters(av_len(max_clusters_av) + 1);
	std::vector<double> node_degree(num_reads + 1, 0.0);
	std::vector<double> sum_internal_edges(num_reads + 1, 0.0);
	std::vector<uint32_t> num_nodes(num_reads + 1, 1);
	std::vector<double> cluster_degree(num_reads + 1, 0.0);
	std::vector<uint32_t> node_queue;
	std::vector<bool> cluster_updated(num_reads + 1, 0);
	std::vector<double> weights;
	std::vector<uint32_t> read_ids;
	uint64_t relayed_seq_graph_index[2];
	double last_Qx = log(0.0);
	
	// 各極大クラスターについて処理
	for (uint32_t i = 0;i < max_clusters.size();i++) {
		SV* max_cluster_sv = *av_fetch(max_clusters_av, i, 0);
		uint32_t* max_cluster = (uint32_t*)SvPV_nolen(max_cluster_sv);
		max_clusters[i].assign(max_cluster, max_cluster + SvCUR(max_cluster_sv) / sizeof(uint32_t));
	}
	
	// 変数を初期化
	for (uint32_t read_id = 1;read_id <= num_reads;read_id++) {
		// シーケンスグラフインデックスが0の場合は以下の処理をスキップ
		if (!seq_graph_index[read_id]) {continue;}
		
		// シーケンスグラフファイルのファイルポインタをシーケンスグラフインデックスの位置にセット
		fseek(NSG, seq_graph_index[read_id], 0);
		
		// シーケンスグラフファイルから16バイト読み込む
		fread(relayed_seq_graph_index, sizeof(uint64_t), 2, NSG);
		
		// リード数に応じて重みリストのサイズを調整
		weights.resize(relayed_seq_graph_index[1] >> 32);
		
		// シーケンスグラフファイルからリード数×8バイト読み込む
		fread(weights.data(), sizeof(double), weights.size(), NSG);
		
		// 重みの和を取得 (自己リンクの1を減算)
		node_degree[read_id] = std::accumulate(weights.begin(), weights.end(), -1.0, [](double x, double y) {return x + fabs(y);});
		cluster_degree[read_id] = node_degree[read_id];
		
		// ノードキューを初期化
		node_queue.push_back(read_id);
	}
	
	// シーケンスグラフインデックスが0でないリード数を取得
	num_reads = node_queue.size();
	
	// グラフ次数を算出
	double graph_degree = std::accumulate(cluster_degree.begin(), cluster_degree.end(), 0.0);
	
	// グラフ密度を算出
	double graph_density = calc_density(graph_degree, num_reads);
	
	// 初期Qxを算出
	double Qx = std::accumulate(cluster_degree.begin(), cluster_degree.end(), 0.0, [graph_degree, graph_density](double x, double y) {return x + calc_partial_Qx(0.0, 1, y, graph_degree, graph_density);});
	
	// メルセンヌ・ツイスターによる擬似乱数生成器を作成
	std::random_device seed;
	std::mt19937 mt(seed());
	
	// Qxが向上した場合は処理を継続
	while (Qx > last_Qx) {
		// 現在のQxを保存
		last_Qx = Qx;
		
		// クラスター更新フラグをリセット
		cluster_updated.assign(cluster_updated.size(), 0);
		
		// 各ノードについて処理
		for (uint32_t i = num_reads - 1;i < num_reads;i--) {
			// 0以上i以下の整数値の一様分布器を作成
			std::uniform_int_distribution<> rnd(0, i);
			
			// ランダムにノードを選択
			uint32_t j = rnd(mt);
			uint32_t node = node_queue[j];
			node_queue[j] = node_queue[i];
			node_queue[i] = node;
			
			// シーケンスグラフファイルのファイルポインタをシーケンスグラフインデックスの位置にセット
			fseek(NSG, seq_graph_index[node], 0);
			
			// シーケンスグラフファイルから16バイト読み込む
			fread(relayed_seq_graph_index, sizeof(uint64_t), 2, NSG);
			
			// リード数に応じてリードIDリスト及び重みリストのサイズを調整
			weights.resize(relayed_seq_graph_index[1] >> 32);
			read_ids.resize(relayed_seq_graph_index[1] >> 32);
			
			// シーケンスグラフファイルからリード数×8バイト読み込む
			fread(weights.data(), sizeof(double), weights.size(), NSG);
			
			// シーケンスグラフファイルのファイルポインタをリレーするシーケンスグラフインデックスの位置にセット
			fseek(NSG, relayed_seq_graph_index[0], 0);
			
			// シーケンスグラフファイルからリード数×4バイト読み込む
			fread(read_ids.data(), sizeof(uint32_t), read_ids.size(), NSG);
			
			// 所属クラスターを取得
			uint32_t cluster = cluster_assignment[node];
			
			// 変数を宣言 (所属クラスターは自己リンクの2をあらかじめ減算)
			std::unordered_map<uint32_t, double> sum_linked_edges{{cluster, -2.0}, {0, 0.0}};
			
			// 隣接クラスターに対するエッジの重みの総和を算出
			for (uint32_t j = 0;j < read_ids.size();j++) {sum_linked_edges[cluster_assignment[read_ids[j]]] += 2.0 * fabs(weights[j]);}
			
			// 所属クラスターから転出する場合のQx変化量を算出
			double basal_Qx_change = calc_partial_Qx(sum_internal_edges[cluster] - sum_linked_edges[cluster], num_nodes[cluster] - 1, cluster_degree[cluster] - node_degree[node], graph_degree, graph_density) - calc_partial_Qx(sum_internal_edges[cluster], num_nodes[cluster], cluster_degree[cluster], graph_degree, graph_density);
			
			// 隣接クラスターに転入 (あるいは所属クラスターから独立) する場合のQx変化量を算出
			std::unordered_map<uint32_t, double> Qx_change;
			for (const auto& sum_linked_edge : sum_linked_edges) {
				Qx_change[sum_linked_edge.first] = calc_partial_Qx(sum_internal_edges[sum_linked_edge.first] + sum_linked_edge.second, num_nodes[sum_linked_edge.first] + 1, cluster_degree[sum_linked_edge.first] + node_degree[node], graph_degree, graph_density) - calc_partial_Qx(sum_internal_edges[sum_linked_edge.first], num_nodes[sum_linked_edge.first], cluster_degree[sum_linked_edge.first], graph_degree, graph_density);
			}
			
			// 所属クラスターを変えない場合のQx変化量を定義
			Qx_change[cluster] = -basal_Qx_change;
			
			// Qx変化量が最大となるクラスターを選出
			uint32_t target_cluster = std::max_element(Qx_change.begin(), Qx_change.end(), [](const auto& x, const auto& y) {return x.second < y.second | x.second == y.second & x.first < y.first;})->first;
			
			// 所属クラスターから独立する場合は新規クラスターを登録
			if (!target_cluster) {
				while (num_nodes[target_cluster]) {target_cluster++;}
				sum_linked_edges[target_cluster] = 0.0;
				Qx_change[target_cluster] = Qx_change[0];
			}
			
			// クラスター内エッジの重みの総和を更新
			sum_internal_edges[target_cluster] += sum_linked_edges[target_cluster];
			sum_internal_edges[cluster] -= sum_linked_edges[cluster];
			
			// ノード数を更新
			num_nodes[target_cluster]++;
			num_nodes[cluster]--;
			
			// クラスター次数を更新
			cluster_degree[target_cluster] += node_degree[node];
			cluster_degree[cluster] -= node_degree[node];
			
			// 所属クラスターを更新
			cluster_assignment[node] = target_cluster;
			
			// Qxを更新
			Qx += Qx_change[target_cluster] + basal_Qx_change;
			
			// クラスター更新フラグをチェック
			cluster_updated[target_cluster] = target_cluster != cluster;
			cluster_updated[cluster] = target_cluster != cluster;
		}
		
		// 各極大クラスターについて属するノードの所属クラスターのいずれかが更新されたか否かで分類
		auto part = std::partition(max_clusters.begin(), max_clusters.end(), [cluster_updated, cluster_assignment](const auto& x) {return std::any_of(x.begin(), x.end(), [=](uint32_t y) {return cluster_updated[cluster_assignment[y]];});});
		
		// 所属クラスターが固定されたノードを取得
		std::unordered_map<uint32_t, uint32_t> fixed_nodes;
		for (auto max_cluster = part;max_cluster != max_clusters.end();max_cluster++) {
			for (uint32_t i = 0;i < (*max_cluster).size();i++) {
				fixed_nodes[(*max_cluster)[i]] = 0;
			}
		}
		
		// 属するノードの所属クラスターがいずれも固定された極大クラスターを削除
		max_clusters.erase(part, max_clusters.end());
		
		// 所属クラスターが固定された各ノードについてノードキューの位置を取得
		for (uint32_t i = 0;i < num_reads;i++) {
			if (fixed_nodes.find(node_queue[i]) != fixed_nodes.end()) {fixed_nodes[node_queue[i]] = i;}
		}
		
		// 固定されたノード数だけリード数を減算
		num_reads -= fixed_nodes.size();
		
		// 変数を宣言
		auto fixed_node = fixed_nodes.begin();
		
		// 所属クラスターが固定された各ノードについてノードキューの末尾に配置
		for (uint32_t i = num_reads;i < num_reads + fixed_nodes.size();i++, fixed_node++) {
			std::swap(node_queue[fixed_node->second], node_queue[i]);
			if (fixed_nodes.find(node_queue[fixed_node->second]) != fixed_nodes.end()) {fixed_nodes[node_queue[fixed_node->second]] = fixed_node->second;}
		}
	}
	
	// NSGファイルを閉じる
	fclose(NSG);
	
	// Qxを返す
	return Qx / graph_degree;
}

__CPP__
// コンセンサスリードシーケンスの構築 main::create_consensus_sequence(リードシーケンスリストリファレンス, 最大アラインメント数, アラインメント方法, マッチ得点, ミスマッチ減点, ギャップ挿入減点, ギャップ伸長減点)
SV* create_consensus_sequence(AV* read_seqs, unsigned int max_aligns, int align_method, int match_award, int mismatch_penalty, int gap_opening_penalty, int gap_extension_penalty) {
	// 塩基コードを定義
	const char* base_code = "ACTGI*-N";
	
	// リード数を取得
	const uint32_t num_reads = av_len(read_seqs);
	
	// 定数を定義
	const __m128i mask = _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF);
	
	// 変数を宣言
	std::vector<std::string> aligned_read_seqs;
	std::string consensus_read_seq;
	
	// アラインメントエンジンを定義
	auto align_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(align_method), match_award, mismatch_penalty, gap_opening_penalty, gap_extension_penalty);
	
	// グラフを作成
	auto graph = spoa::createGraph();
	
	// 各リードシーケンスについてアラインメントを実行
	for (int i = 0;i <= num_reads && i <= max_aligns;i++) {
		auto read_seq = SvPV_nolen(*av_fetch(read_seqs, i, 0));
		auto alignment = align_engine->align(read_seq, graph);
		graph->add_alignment(alignment, read_seq);
	}
	
	// マルチプルアラインメントを作成
	graph->generate_multiple_sequence_alignment(aligned_read_seqs);
	
	// アラインメント長を取得
	const uint32_t align_len = aligned_read_seqs[0].size();
	
	// 変数を宣言
	std::vector<uint16_t> depth(align_len * 8);
	
	// 各アラインメント済みリードシーケンスについて処理
	for (const auto& aligned_seq : aligned_read_seqs) {
		// 各サイトについて該当する塩基の深度を集計
		for (uint32_t i = 0;i < align_len;i++) {depth[i * 8 + (aligned_seq.c_str()[i] >> 1 & 0x07)]++;}
	}
	
	// 各サイトについて処理
	for (uint32_t i = 0;i < align_len;i++) {
		// 多数決で塩基を決定
		consensus_read_seq += base_code[_mm_cvtsi128_si32(_mm_minpos_epu16(_mm_xor_si128(_mm_load_si128(reinterpret_cast<__m128i*>(&depth[i * 8])), mask))) >> 16];
	}
	
	// コンセンサスシーケンスを返す
	return newSVpvf("%s", consensus_read_seq.c_str());
}
