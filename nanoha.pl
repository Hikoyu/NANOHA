#!/usr/bin/env perl

# Copyright (c) 2020 Hikoyu Suzuki
# This software is released under the MIT License.

use strict;
use warnings;
use Getopt::Std;
use threads;

# ソフトウェアを定義
### 編集範囲 開始 ###
my $software = "nanoha.pl";	# ソフトウェアの名前
my $version = "ver.1.3.0";	# ソフトウェアのバージョン
my $note = "NANOHA is Network-based Assortment of Noisy On-target reads for High-accuracy Alignments.\n  This software assorts on-target PacBio/Nanopore reads such as target amplicon sequences.";	# ソフトウェアの説明
my $usage = "<required items> [optional items]";	# ソフトウェアの使用法 (コマンド非使用ソフトウェアの時に有効)
### 編集範囲 終了 ###

# コマンドを定義
my %command;
### 編集範囲 開始 ###
$command{"sketch"} = "Sketch out sequence reads";
$command{"assort"} = "Assort sequence reads";
$command{"unify"} = "\tUnify sequence reads in the same cluster";
$command{"convert"} = "Convert sequence reads from FASTA format to FASTQ format";
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
use List::MoreUtils 0.31_01;
use Math::BigFloat;
use Number::AnyBase 1.60000;
use Inline (C => Config => CC => exists($ENV{"CC"}) ? $ENV{"CC"} : 'cc');
use Inline (C => 'DATA', NAME => 'NANOHA::LLCS::NOVEC', CCFLAGS => '-std=c99');
use Inline (C => Config => AUTO_INCLUDE => '#include <x86intrin.h>');
use Inline (C => 'DATA', NAME => 'NANOHA::LLCS::SSE41', CCFLAGS => '-std=c99 -msse4.1');
use Inline (C => 'DATA', NAME => 'NANOHA::LLCS::AVX2', CCFLAGS => '-std=c99 -mavx2');
use Inline (C => 'DATA', NAME => 'NANOHA::LLCS::AVX512', CCFLAGS => '-std=c99 -mavx512f -mavx512bw');
use Inline (CPP => 'DATA', NAME => 'NANOHA::MSA', AUTO_INCLUDE => ['#include <x86intrin.h>', '#include <spoa/spoa.hpp>'], CC => exists($ENV{"CXX"}) ? $ENV{"CXX"} : 'g++', CCFLAGS => '-std=c++11 -march=native', INC => exists($ENV{"SPOA_INC"}) ? "-I$ENV{SPOA_INC}" : '-I/usr/local/inculde', LIBS => [exists($ENV{"SPOA_LIB"}) ? "-Wl,-rpath,$ENV{SPOA_LIB} -L$ENV{SPOA_LIB}" : '-Wl,-rpath,/usr/local/lib -L/usr/local/lib', '-lspoa']);
no warnings 'portable';

# 定数を定義
use constant {
	max_num_reads => 4294967295,	# max_num_reads => リード数上限値
	max_depth => 65535,				# max_depth => 深度上限値
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
	# オプションを追加
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
	
	# 入力ファイルを確認
	&exception::error("input file not specified") unless @ARGV or -p STDIN;
	&common::check_files(\@ARGV);
	
	# 共有変数を定義
	my $seq_index : shared;
	
	# プロセスIDとプログラム開始時刻をファイルヘッダーとしてシーケンスインデックスに登録
	$seq_index = pack("NN", 0, 0);
	
	# Nanoha Sequence Read (NSR) ファイルを作成
	open(NSR, ">", "$opt{o}.nsr") or &exception::error("failed to make file: $opt{o}.nsr");
	
	# ファイルヘッダーをバイナリ形式でNSRファイルに出力
	print NSR $seq_index;
	
	# 変数を宣言
	my @worker_thread = ();
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
		$worker_thread[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューを読み込みながら処理
			while (defined(my $dat = $input->dequeue)) {
				# 変数を宣言
				my %kmers = ();
				
				# リードIDとリードシーケンスをデシリアライズして取得
				my ($read_id, $read_seq) = unpack("NA*", $dat);
				
				# 順鎖のシーケンスに含まれるk-merを取得
				$kmers{$_}->[0]++ foreach map {substr($read_seq, $_, $opt{"k"})} 0..length($read_seq) - $opt{"k"};
				
				# 逆鎖のシーケンスに含まれるk-merを取得 (-s未指定時)
				$kmers{$_}->[0]++ foreach map {&common::complementary($_)} $opt{"s"} ? () : keys(%kmers);
				
				# 基数変換したk-merをシード値に用いて指定した個数のハッシュ値を生成
				srand($convert->to_dec($_)) and push(@{$kmers{$_}}, map {rand} 1..$opt{"n"}) foreach keys(%kmers);
				
				# minimizerを生成
				my @minimizers = map {List::Util::reduce {$kmers{$a}->[$_] < $kmers{$b}->[$_] ? $a : $b} keys(%kmers)} 1..$opt{"n"};
				
				# リードID及び各minimizerとそのカウント数をシリアライズして出力キューに追加
				$output->enqueue(pack("N(Nn)*", $read_id, map {$convert->to_dec($_), $kmers{$_}->[0]} @minimizers));
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# データストリームスレッドを作成
	## ここからデータストリームスレッドの処理 ##
	$worker_thread[$opt{"p"}] = threads::async {
		# このスレッドのみを終了可能に変更
		threads->set_thread_exit_only(1);
		
		# Nanoha Sequence Sketch (NSS) ファイルを作成
		open(NSS, ">", "$opt{o}.nss") or &exception::error("failed to make file: $opt{o}.nss");
		
		# ファイルヘッダーをバイナリ形式でNSSファイルに出力
		print NSS $seq_index;
		
		# k-merサイズ及びストランド特異性フラグとminimizerの生成個数をバイナリ形式でNSSファイルに出力
		print NSS pack("CC", !!$opt{"s"} << 4 | $opt{"k"}, $opt{"n"});
		
		# 変数を宣言
		my $seq_sketch_index = "";
		
		# 出力キューを読み込みながら処理
		while (defined(my $dat = $output->dequeue)) {
			# NSSファイルのファイルポインタの位置をシーケンススケッチインデックスに登録
			vec($seq_sketch_index, vec($dat, 0, 32) * 2, 64) = tell(NSS);
			
			# 各minimizerとそのカウント数をバイナリ形式でNSSファイルに出力
			print NSS substr($dat, 4);
		}
		
		# NSSファイルを閉じる
		close(NSS);
		
		# シーケンススケッチインデックスをシーケンスインデックスに統合
		$seq_index |= $seq_sketch_index;
		
		# スレッドを終了
		return(1);
	};
	## ここまでデータストリームスレッドの処理 ##
	
	# 変数を宣言
	my $read_id = 0;
	my $seq_read_index = "";
	
	# NSRファイルのファイルポインタの位置をシーケンスリードインデックスに登録
	vec($seq_read_index, 1, 64) = tell(NSR);
	
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
		
		# 実行中のスレッド数が指定値+1と一致している場合はリードIDとリードシーケンスをシリアライズして入力キューに追加
		$input->enqueue(pack("NA*", $read_id, $read_seq)) if threads->list(threads::running) == $opt{"p"} + 1;
		
		# リードシーケンスを2進数に変換
		$read_seq =~ s/A/00/g;
		$read_seq =~ s/C/10/g;
		$read_seq =~ s/G/01/g;
		$read_seq =~ s/T/11/g;
		
		# リード長とリードシーケンスをバイナリ形式でNSRファイルに出力
		print NSR pack("Nb*", $read_len, $read_seq);
		
		# NSRファイルのファイルポインタの位置をシーケンスインデックスに登録
		vec($seq_read_index, $read_id * 2 + 1, 64) = tell(NSR);
		
		# リードIDが上限値に達した場合、読み込みを終了
		last if $read_id == $opt{"u"};
	}
	
	# NSRファイルを閉じる
	close(NSR);
	
	# シーケンスリードインデックスをシーケンスインデックスに統合
	$seq_index |= $seq_read_index;
	
	# 入力キューを終了
	$input->end;
	
	# 並列処理の各ワーカースレッドが終了するまで待機
	$worker_thread[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 出力キューを終了
	$output->end;
	
	# データストリームスレッドが終了するまで待機
	$worker_thread[$opt{"p"}]->join or &exception::error("data stream thread abnormally exited");
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	# Nanoha Sequence Index (NSI) ファイルを作成
	open(NSI, ">", "$opt{o}.nsi") or &exception::error("failed to make file: $opt{o}.nsi");
	
	# ファイルヘッダー及びシーケンスインデックスをバイナリ形式でNSIファイルに出力
	print NSI $seq_index;
	
	# NSIファイルを閉じる
	close(NSI);
	
	return(1);
}

# サブルーチンを追加

## ここからassortコマンドのパッケージ ##
package assort;

# コマンドとオプションを定義
sub define {
	$note = "Assort sequence reads under specified conditions.";
	$usage = "<prefix>";
	$option{"p INT "} = "Number of parallel worker threads <1-> [1]";
	$option{"t INT "} = "Number of trials per worker thread for clustering sequence reads <1-> [1]";
	$option{"m INT "} = "Maximum number of k-mer minimizers to be used from each sequence read <1-255> [255]";
	$option{"w INT "} = "Word size index for calculating LLCS in Perl XS code <0-3>";
	# オプションを追加
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT >= 1: -t $opt{t}") if $opt{"t"} !~ /^\d+$/ or $opt{"t"} < 1;
	&exception::error("specify INT 1-255: -m $opt{m}") if $opt{"m"} !~ /^\d+$/ or $opt{"m"} < 1 or $opt{"m"} > 255;
	&exception::error("specify INT 0-3: -w $opt{w}") if defined($opt{"w"}) and ($opt{"w"} !~ /^\d+$/ or $opt{"w"} < 0 or $opt{"w"} > 3);
	
	# LLCSの計算エンジンを表示
	print STDERR "Use ", defined($opt{"w"}) ? "Perl XS code" : "Pure Perl code", " for LLCS calculation\n";
	
	# 入力ファイルのプレフィックスを取得
	my $input_prefix = shift(@ARGV);
	
	# 入力ファイルのプレフィックスを確認
	&exception::error("input file prefix not specified") unless defined($input_prefix);
	
	# 入力ファイルを確認
	&common::check_files(["$input_prefix.nsi", "$input_prefix.nss", "$input_prefix.nsr"]);
	
	# Nanoha Sequence Index (NSI) ファイルを開く
	open(NSI, "<", "$input_prefix.nsi") or &exception::error("failed to open file: $input_prefix.nsi");
	
	# NSIファイルを読み込む
	print STDERR "Loading sequence indexes...";
	read(NSI, my $seq_index, -s "$input_prefix.nsi");
	print STDERR "completed\n";
	
	# NSIファイルを閉じる
	close(NSI);
	
	# Nanoha Sequence Sketch (NSS) ファイルを開く
	open(NSS, "<", "$input_prefix.nss") or &exception::error("failed to open file: $input_prefix.nss");
	
	# NSSファイルを先頭から10バイト読み込む
	read(NSS, my $info, 10);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsi, $input_prefix.nss") unless vec($seq_index, 0, 64) == vec($info, 0, 64);
	
	# ストランド特異性フラグを取得
	my $strand_speicificity = vec($info, 17, 4);
	
	# minimizerの使用個数を定義
	my $num_minimizers = List::Util::min(vec($info, 9, 8), $opt{"m"});
	
	# 変数を宣言
	my @minimizer_groups = ();
	my $read_id = 0;
	
	# 各リードについてシーケンススケッチインデックスを取得しながら処理
	print STDERR "Loading sequence sketches...";
	while (my $seq_sketch_index = vec($seq_index, ($read_id + 1) * 2, 64)) {
		# 変数を宣言
		my $minimizer = "";
		
		# リードIDを更新
		$read_id++;
		
		# NSSファイルのファイルポインタをシーケンススケッチインデックスの位置にセット
		seek(NSS, $seq_sketch_index, 0);
		
		# minimizerの使用個数分だけNSSファイルから6バイト読み込みminimizerグループにリードIDをシリアライズして追加
		read(NSS, $minimizer, 6) and $minimizer_groups[$_]->{$minimizer} .= pack("N", $read_id) foreach 0..$num_minimizers - 1;
	}
	print STDERR "completed\n";
	
	# NSSファイルを閉じる
	close(NSS);
	
	### MinHash法に基づくリード間のJaccard類似度の算出 ###
	# 共有変数を定義
	my @edge_list : shared;
	@edge_list = ("") x ($read_id + 1);
	
	# 変数を宣言
	my @worker_thread = ();
	my $num_error_threads = 0;
	
	# 入力キューを作成
	my $queue = Thread::Queue->new;
	
	# 入力キューの要素数上限を定義
	$queue->limit = $opt{"p"};
	
	# 指定されたワーカースレッド数で並列処理
	print STDERR "Calculating Jaccard similarities based on MinHash...";
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_thread[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューを読み込みながら処理
			while (defined(my $dat = $queue->dequeue)) {
				# minimizerグループをデシリアライズ
				my @minimizer_group = unpack("N*", $dat);
				
				# minimizerグループ内の各リードについて処理
				while (my $read_id = shift(@minimizer_group)) {
					# エッジをデシリアライズして取得
					my %edges = unpack("(NC)*", $edge_list[$read_id]);
					
					# エッジの重みを加算
					$edges{$_}++ foreach @minimizer_group;
					
					# シリアライズしてエッジリストに登録
					$edge_list[$read_id] = pack("(NC)*", %edges);
				}
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# minimizerの使用個数分だけ処理
	while (my $minimizer_group = shift(@minimizer_groups)) {
		# 各minimizerグループについて実行中のスレッド数が指定値と一致している場合はグループ内のリードIDをシリアライズしたまま入力キューに追加
		threads->list(threads::running) == $opt{"p"} and $queue->enqueue($_) foreach values(%{$minimizer_group});
	}
	
	# 入力キューを終了
	$queue->end;
	
	# 並列処理の各ワーカースレッドが終了するまで待機
	$worker_thread[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	### LLCSに基づくリード間のJaccard類似度の算出 ###
	# 追加ヌル文字数を算出
	my $additional_nulls = defined($opt{"w"}) ? (0x10 << $opt{"w"}) - 1 : 0x0F;
	
	# 変数を宣言
	# 入出力キューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 入力キューの要素数上限を定義
	$input->limit = $opt{"p"};
	
	# 指定されたワーカースレッド数で並列処理
	print STDERR "Calculating Jaccard similarities based on LLCS...";
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_thread[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# Nanoha Sequence Read (NSR) ファイルを開く
			open(NSR, "<", "$input_prefix.nsr") or &exception::error("failed to open file: $input_prefix.nsr");
			
			# NSRファイルを先頭から8バイト読み込む
			read(NSR, my $file_header, 8);
			
			# ファイルヘッダーの一致を確認
			&exception::error("file headers not matched: $input_prefix.nsi, $input_prefix.nsr") unless vec($seq_index, 0, 64) == vec($file_header, 0, 64);
			
			# 入力キューを読み込みながら処理
			while (defined(my $dat = $input->dequeue)) {
				# 変数を宣言
				my @read_seqs = ();
				my @read_lens = ();
				
				# リードIDをデシリアライズして取得
				my @read_ids = unpack("N*", $dat);
				
				# 各リードについて処理
				foreach my $read_id (@read_ids) {
					# シーケンスリードインデックスを取得
					my $seq_read_index = vec($seq_index, ($read_id - 1) * 2 + 1, 64);
					
					# NSRファイルのファイルポインタをシーケンスリードインデックスの位置にセット
					seek(NSR, $seq_read_index, 0);
					
					# NSRファイルから4バイト読み込みリード長リストに追加
					read(NSR, my $read_len, 4);
					
					# NSRファイルからリードシーケンスをバイナリ形式で読み込む
					read(NSR, my $read_seq, vec($seq_index, $read_id * 2 + 1, 64) - $seq_read_index - 4);
					
					# リードシーケンスリストにリードシーケンスをバイナリ形式で追加
					push(@read_seqs, $read_seq . chr(0) x $additional_nulls);
					
					# リード長リストにリード長を追加
					push(@read_lens, vec($read_len, 0, 32));
				}
				
				# LLCSに基づくJaccard類似度を算出
				my @LLCS_Jaccard_similarities = map {$_ * 127} &calc_LLCS_Jaccard_similarities(\@read_seqs, \@read_lens, $strand_speicificity, $opt{"w"});
				
				# リードIDとLLCSに基づくJaccard類似度をシリアライズして出力キューに追加
				$output->enqueue(pack("N(Nc)*", shift(@read_ids), List::MoreUtils::mesh(@read_ids, @LLCS_Jaccard_similarities)));
			}
			
			# NSRファイルを閉じる
			close(NSR);
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 各リードについて実行中のスレッド数が指定値と一致している場合は自身のリードID及びエッジリストに登録されているリードIDをシリアライズして入力キューに追加
	threads->list(threads::running) == $opt{"p"} and $input->enqueue(pack("N*", $_, unpack("(Nx)*", $edge_list[$_]))) foreach 1..$read_id;
	
	# 入力キューを終了
	$input->end;
	
	# 並列処理の各ワーカースレッドが終了するまで待機
	$worker_thread[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 出力キューを終了
	$output->end;
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	# 変数を宣言
	my @LLCS_Jaccard_similarities = ();
	
	# 出力キューからLLCSに基づくJaccard類似度を読み込み登録
	print STDERR "Modifying the weight of each edge...";
	while (defined(my $dat = $output->dequeue)) {$LLCS_Jaccard_similarities[vec($dat, 0, 32)] = substr($dat, 4);}
	
	# 各リードについて処理
	for (my $i = 1;$i <= $read_id;$i++) {
		# LLCSに基づくJaccard類似度をデシリアライズ
		my %LLCS_Jaccard_similarity = unpack("(Nc)*", $LLCS_Jaccard_similarities[$i]);
		
		# エッジの重みを修正して登録
		$edge_list[$i] = pack("(Nc)*", List::Util::pairmap {$a => $LLCS_Jaccard_similarity{$a} * $b / $num_minimizers} unpack("(NC)*", $edge_list[$i]));
	}
	
	# ペアを入れ替えてエッジを登録
	map {$edge_list[$_->[0]] .= $_->[1]} List::Util::pairmap {[$a, pack("Nc", $_, $b)]} List::Util::pairgrep {$a > $_} unpack("(Nc)*", $edge_list[$_]) foreach 1..$read_id;
	print STDERR "completed\n";
	
	### excess modularity density (Qx) に基づくリードのクラスタリング ###
	# 共有変数を定義
	my $max_Qx : shared;
	my $cluster_assignment : shared;
	
	# Qx最大値の初期値を設定
	$max_Qx = -Inf;
	
	# 指定されたワーカースレッド数で並列処理
	print STDERR "Creating sequence clusters...";
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_thread[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 指定した回数に達するまで処理
			for (my $num_trials = 0;$num_trials < $opt{"t"};$num_trials++) {
				# excess modularity density (Qx) に基づいてリードをクラスタリング
				my ($Qx, @cluster_assignment) = &create_sequence_clusters(\@edge_list);
				
				# 共有変数をロック
				lock($max_Qx);
				lock($cluster_assignment);
				
				# Qxが最大値を超えた場合はQx最大値及びクラスター割り当てリストを更新してトライアル回数をリセット
				($max_Qx, $cluster_assignment, $num_trials) = ($Qx, pack("N*", @cluster_assignment), -1) if $Qx > $max_Qx;
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 並列処理の各ワーカースレッドが終了するまで待機
	$worker_thread[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	### シーケンス方向の決定 ###
	# 変数を宣言
	my %clusters = ();
	my $seq_ori = substr($seq_index, 0, 8);
	
	# 各リードについて所属するクラスターに登録
	push(@{$clusters{vec($cluster_assignment, $_, 32)}}, $_) foreach 1..$read_id;
	
	# クラスターサイズリストを作成
	my @cluster_sizes = map {scalar(@{$_})} values(%clusters);
	
	# 判別分析法（大津の二値化）に基づく深度の推奨閾値を算出
	my $cutoff_depth = List::Util::reduce {$a->[1] > $b->[1] ? $a : $b}
	List::Util::pairmap {[List::Util::min(@{$b}), defined($a) ? @{$a} * @{$b} * (List::Util::sum(@{$a}) / @{$a} - List::Util::sum(@{$b}) / @{$b}) ** 2 : 0]}
	map {my $c = $_;List::MoreUtils::part {$_ >= $c} @cluster_sizes}
	sort {$a <=> $b} List::Util::uniq(@cluster_sizes);
	
	# ストランド特異性フラグを登録
	vec($seq_ori, 64, 1) = $strand_speicificity;
	
	# Nanoha Sequence Cluster (NSC) ファイルを作成
	open(NSC, ">", "$input_prefix.nsc") or &exception::error("failed to make file: $input_prefix.nsc");
	
	# ファイルヘッダーをバイナリ形式でNSCファイルに出力
	print NSC substr($seq_index, 0, 8);
	
	# 深度の推奨閾値をバイナリ形式でNSCファイルに出力
	print NSC pack("N", $cutoff_depth->[0]);
	
	# 各クラスターについて処理
	print STDERR "Determining sequence orientations...";
	foreach my $cluster (sort {$a <=> $b} keys(%clusters)) {
		# クラスター密度と所属リードのリードIDをバイナリ形式でNSCファイルに出力
		print NSC pack("N/N*", @{$clusters{$cluster}});
		
		# 先頭リードのリードIDを取得してリードキューに追加
		my @read_queue = shift(@{$clusters{$cluster}});
		
		# 残りの所属リードのリードIDをメンバーに登録
		my %members = map {$_ => 1} @{$clusters{$cluster}};
		
		# リードキューからリードIDを取得して処理
		while (my $read_id = shift(@read_queue)) {
			# 変数を宣言
			my %edges = unpack("(Nc)*", $edge_list[$read_id]);
			
			# エッジリストに登録されている各リードのうちメンバーに登録されているものについてメンバーから削除するとともにリードキューに追加してシーケンス方向を決定
			delete($members{$_}) and push(@read_queue, $_) and vec($seq_ori, $_ + 64, 1) = vec($seq_ori, $read_id + 64, 1) ^ ($edges{$_} < 0) foreach keys(%edges);
		}
	}
	print STDERR "completed\n";
	
	# NSCファイルを閉じる
	close(NSC);
	
	# Nanoha Sequence Orientation (NSO) ファイルを作成
	open(NSO, ">", "$input_prefix.nso") or &exception::error("failed to make file: $input_prefix.nso");
	
	# ファイルヘッダー及びシーケンス方向をバイナリ形式でNSOファイルに出力
	print NSO $seq_ori;
	
	# NSOファイルを閉じる
	close(NSO);
	
	# クラスター数を表示
	print STDERR scalar(@cluster_sizes), " sequence cluster", @cluster_sizes > 1 ? "s" : "", " found\n";
	
	# Qx最大値を表示
	print STDERR "Excess modularity density: $max_Qx\n";
	
	# 深度の推奨閾値を表示
	print STDERR "Recommended cutoff depth: $cutoff_depth->[0]\n" if defined($cutoff_depth);
	
	return(1);
}

# LLCSに基づくJaccard類似度を算出 classify::calc_LLCS_Jaccard_similarities(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性フラグ, ワードサイズ指数)
sub calc_LLCS_Jaccard_similarities {
	# 引数を取得
	my ($read_seqs, $read_lens, $strand_speicificity, $word_size_index) = @_;
	
	# ワードサイズ指数が定義されている場合はSIMD演算を用いたCの関数でLLCSに基づくJaccard類似度を算出して返す
	return(&{\&{"main::calc_LLCS_Jaccard_similarities_${word_size_index}"}}($read_seqs, $read_lens, $strand_speicificity)) if defined($word_size_index);
	
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
		my @v = ();
		my $u = 0;
		my $w = 0;
		
		# 計算要素を初期化
		@v = (0xFFFFFFFFFFFFFFFF) x $num_query_blocks;
		push(@v, 0xFFFFFFFFFFFFFFFF) if @v % 4;
		push(@v, 0xFFFFFFFFFFFFFFFF) if @v % 4;
		push(@v, 0xFFFFFFFFFFFFFFFF) if @v % 4;
		
		# 順鎖に対するLLCSを算出
		for (my $i = 0;$i < $num_blocks;$i++) {
			my $end = $i + !!($read_len & 0x1F) < $num_blocks ? 64 : $read_len << 1 & 0x3E;
			for (my $j = 0;$j < $end;$j += 2) {
				my $base = $read_seq->[$i] >> $j & 0x03;
				my $carry = 0;
				for (my $k = 0;$k < $num_query_blocks;$k++) {
					$u = $v[$k] & $query_matrix[$base]->[$k];
					$w = ($v[$k] & 0x7FFFFFFFFFFFFFFF) + ($u & 0x7FFFFFFFFFFFFFFF) + $carry;
					$w ^= ($v[$k] ^ $u) & 0x8000000000000000;
					$carry = $w < $v[$k];
					$v[$k] = $w | ($v[$k] - $u);
				}
			}
		}
		my $for_llcs = List::Util::sum(
			map {($_ & 0x00000000FFFFFFFF) + ($_ >> 32 & 0x00000000FFFFFFFF)}
			map {($_ & 0x0000FFFF0000FFFF) + ($_ >> 16 & 0x0000FFFF0000FFFF)}
			List::Util::pairmap {$a + $b}
			map {($_ & 0x00FF00FF00FF00FF) + ($_ >> 8 & 0x00FF00FF00FF00FF)}
			map {($_ & 0x0F0F0F0F0F0F0F0F) + ($_ >> 4 & 0x0F0F0F0F0F0F0F0F)}
			List::Util::pairmap {$a + $b}
			map {($_ & 0x3333333333333333) + ($_ >> 2 & 0x3333333333333333)}
			map {($_ & 0x5555555555555555) + ($_ >> 1 & 0x5555555555555555)}
			map {~$_} @v
		);
		
		# ストランド特異性フラグが立っている場合は順鎖に対するLLCSに基づくJaccard類似度を算出して以下の処理をスキップ
		push(@LLCS_Jaccard_similarities, $for_llcs / ($query_len + $read_len - $for_llcs)) and next if $strand_speicificity;
		
		# リードシーケンスを相補的な塩基に変換
		$_ = ~$_ foreach @{$read_seq};
		
		# 計算要素を初期化
		@v = (0xFFFFFFFFFFFFFFFF) x $num_query_blocks;
		push(@v, 0xFFFFFFFFFFFFFFFF) if @v % 4;
		push(@v, 0xFFFFFFFFFFFFFFFF) if @v % 4;
		push(@v, 0xFFFFFFFFFFFFFFFF) if @v % 4;
		
		# 逆鎖に対するLLCSを算出
		for (my $i = $num_blocks - 1;$i >= 0;$i--) {
			my $end = $i + !!($read_len & 0x1F) < $num_blocks ? 64 : $read_len << 1 & 0x3E;
			for (my $j = $end - 2;$j >= 0;$j -= 2) {
				my $base = $read_seq->[$i] >> $j & 0x03;
				my $carry = 0;
				for (my $k = 0;$k < $num_query_blocks;$k++) {
					$u = $v[$k] & $query_matrix[$base]->[$k];
					$w = ($v[$k] & 0x7FFFFFFFFFFFFFFF) + ($u & 0x7FFFFFFFFFFFFFFF) + $carry;
					$w ^= ($v[$k] ^ $u) & 0x8000000000000000;
					$carry = $w < $v[$k];
					$v[$k] = $w | ($v[$k] - $u);
				}
			}
		}
		my $rev_llcs = List::Util::sum(
			map {($_ & 0x00000000FFFFFFFF) + ($_ >> 32 & 0x00000000FFFFFFFF)}
			map {($_ & 0x0000FFFF0000FFFF) + ($_ >> 16 & 0x0000FFFF0000FFFF)}
			List::Util::pairmap {$a + $b}
			map {($_ & 0x00FF00FF00FF00FF) + ($_ >> 8 & 0x00FF00FF00FF00FF)}
			map {($_ & 0x0F0F0F0F0F0F0F0F) + ($_ >> 4 & 0x0F0F0F0F0F0F0F0F)}
			List::Util::pairmap {$a + $b}
			map {($_ & 0x3333333333333333) + ($_ >> 2 & 0x3333333333333333)}
			map {($_ & 0x5555555555555555) + ($_ >> 1 & 0x5555555555555555)}
			map {~$_} @v
		);
		
		# LLCSに基づくJaccard類似度を算出 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		push(@LLCS_Jaccard_similarities, $for_llcs >= $rev_llcs ? $for_llcs / ($query_len + $read_len - $for_llcs) : -$rev_llcs / ($query_len + $read_len - $rev_llcs));
	}
	
	# LLCSに基づくJaccard類似度を返す
	return(@LLCS_Jaccard_similarities);
}

# 密度を算出 classify::calc_density(内部エッジの重みの総和, ノード数)
sub calc_density {
	# 引数を取得
	my ($sum_internal_edges, $num_nodes) = @_;
	
	# 密度を算出
	my $density = $num_nodes > 1 ? $sum_internal_edges / $num_nodes / ($num_nodes - 1) : 0;
	
	# 密度を返す
	return($density);
}

# 部分的なQxを算出 classify::calc_partial_Qx(内部エッジの重みの総和, ノード数, クラスター次数, グラフ次数, グラフ密度)
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

# excess modularity density (Qx) に基づくリードのクラスタリング classify::create_sequence_clusters(エッジリストリファレンス)
sub create_sequence_clusters {
	# 引数を取得
	my ($edge_list) = @_;
	
	# リード数を取得
	my $num_reads = @{$edge_list} - 1;
	
	# 変数を宣言
	my @node_degree = map {List::Util::sum0(map {abs($_ / 127)} unpack("(x4c)*", $_))} @{$edge_list};
	my @sum_internal_edges = (0) x @{$edge_list};
	my @num_nodes = (1) x @{$edge_list};
	my @cluster_degree = @node_degree;
	my @cluster_assignment = 0..$num_reads;
	my @missing_clusters = ();
	my $last_Qx = -Inf;
	
	# グラフ次数を算出
	my $graph_degree = List::Util::sum(@node_degree);
	
	# グラフ密度を算出
	my $graph_density = &calc_density($graph_degree, $num_reads);
	
	# 初期Qxを算出
	my $Qx = List::Util::sum(map {&calc_partial_Qx($sum_internal_edges[$_], $num_nodes[$_], $cluster_degree[$_], $graph_degree, $graph_density)} 1..$num_reads);
	
	# Qxが向上した場合は処理を継続
	while ($Qx > $last_Qx) {
		# 現在のQxを保存
		$last_Qx = $Qx;
		
		# ランダムな順序でノードキューを作成
		my @node_queue = List::Util::shuffle(1..$num_reads);
		
		# 各ノードについて処理
		foreach my $node (@node_queue) {
			# エッジリストをデシリアライズ
			my %edges = unpack("(Nc)*", $edge_list->[$node]);
			
			# 所属クラスターを取得
			my $cluster = $cluster_assignment[$node];
			
			# 変数を宣言
			my %sum_linked_edges = ($cluster => 0);
			
			# 新規クラスターを取得
			my $new_cluster = shift(@missing_clusters);
			
			# 新規クラスターが得られた場合
			$sum_linked_edges{$new_cluster} = 0 if defined($new_cluster);
			
			# 隣接クラスターに対するエッジの重みの総和を算出
			$sum_linked_edges{$cluster_assignment[$_]} += 2 * abs($edges{$_} / 127) foreach keys(%edges);
			
			# 現在クラスターから転出した場合のQx変化量を算出
			my $basal_Qx_change = &calc_partial_Qx($sum_internal_edges[$cluster] - $sum_linked_edges{$cluster}, $num_nodes[$cluster] - 1, $cluster_degree[$cluster] - $node_degree[$node], $graph_degree, $graph_density) - &calc_partial_Qx($sum_internal_edges[$cluster], $num_nodes[$cluster], $cluster_degree[$cluster], $graph_degree, $graph_density);
			
			# 隣接クラスターに転入した場合のQx変化量を算出
			my %Qx_change = map {$_ => &calc_partial_Qx($sum_internal_edges[$_] + $sum_linked_edges{$_}, $num_nodes[$_] + 1, $cluster_degree[$_] + $node_degree[$node], $graph_degree, $graph_density) - &calc_partial_Qx($sum_internal_edges[$_], $num_nodes[$_], $cluster_degree[$_], $graph_degree, $graph_density)} keys(%sum_linked_edges);
			
			# Qx変化量が最大となるクラスターを選出
			my $target_cluster = List::Util::reduce {$Qx_change{$a} > $Qx_change{$b} ? $a : $b} sort {$a <=> $b} keys(%Qx_change);
			
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
			
			# 欠番クラスターを追加
			push(@missing_clusters, $new_cluster) if defined($new_cluster) and $new_cluster != $target_cluster;
			push(@missing_clusters, $cluster) unless $num_nodes[$cluster];
			
			# Qxを更新
			$Qx += $Qx_change{$target_cluster} + $basal_Qx_change;
		}
	}
	
	# Qx及びクラスター割り当てリストを返す
	return($Qx / $graph_degree, @cluster_assignment);
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
	&exception::error("specify INT >= 1: -d $opt{d}") if defined($opt{"d"}) and ($opt{"d"} !~ /^\d+$/ or $opt{"d"} < 1);
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT >= 1: -g $opt{g}") if $opt{"g"} !~ /^\d+$/ or $opt{"g"} < 1;
	&exception::error("specify INT <= 1: -e $opt{e}") if $opt{"e"} !~ /^\d+$/ or $opt{"e"} < 1;
	&exception::error("specify INT >= 1: -m $opt{m}") if $opt{"m"} !~ /^\d+$/ or $opt{"m"} < 1;
	&exception::error("specify INT >= 0: -n $opt{n}") if $opt{"n"} !~ /^\d+$/;
	&exception::error("specify INT 1-" . main::max_depth . ": -a $opt{a}") if $opt{"a"} !~ /^\d+$/ or $opt{"a"} < 1 or $opt{"a"} > main::max_depth;
	&exception::error("specify FLOAT 0-1: -q $opt{q}") if $opt{"q"} !~ /^\d+$|^\d+\.\d+$|^\d+[eE]-?\d+$|^\d+\.\d+[eE]-?\d+$/ or $opt{"q"} > 1;
	$opt{"c"} = lc($opt{"c"}) if $opt{"c"};
	
	# 入力ファイルのプレフィックスを取得
	my $input_prefix = shift(@ARGV);
	
	# 入力ファイルのプレフィックスを確認
	&exception::error("input file prefix not specified") unless defined($input_prefix);
	
	# 入力ファイルを確認
	&common::check_files(["$input_prefix.nsc", "$input_prefix.nsi", "$input_prefix.nso", "$input_prefix.nsr"]);
	
	# アラインメント方法の対応番号を定義
	my %align_method = ("local" => 0, "global" => 1, "semi-global" => 2);
	
	# アラインメントに用いるギャップモデルを表示 (-c指定時)
	print STDERR "Use ", $opt{"g"} > $opt{"e"} ? "affine" : "linear", " gap model for the multiple alignment\n" if $opt{"c"};
	
	# Nanoha Sequence Cluster (NSC) ファイルを開く
	open(NSC, "<", "$input_prefix.nsc") or &exception::error("failed to open file: $input_prefix.nsc");
	
	# NSCファイルを先頭から8バイト読み込む
	read(NSC, my $file_header, 8);
	
	# NSCファイルを4バイト読み込む
	read(NSC, my $cutoff_depth, 4);
	
	# 深度の閾値が未定義の場合はNSCファイルから読み込んだ値を使用
	$opt{"d"} = vec($cutoff_depth, 0, 32) and print STDERR "Cutoff depth: $opt{d}\n" unless defined($opt{"d"});
	
	# Nanoha Sequence Index (NSI) ファイルを開く
	open(NSI, "<", "$input_prefix.nsi") or &exception::error("failed to open file: $input_prefix.nsi");
	
	# NSIファイルを読み込む
	print STDERR "Loading sequence indexes...";
	read(NSI, my $seq_index, -s "$input_prefix.nsi");
	print STDERR "completed\n";
	
	# NSIファイルを閉じる
	close(NSI);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsi, $input_prefix.nsc") unless vec($seq_index, 0, 64) == vec($file_header, 0, 64);
	
	# Nanoha Sequence Orientation (NSO) ファイルを開く
	open(NSO, "<", "$input_prefix.nso") or &exception::error("failed to open file: $input_prefix.nso");
	
	# NSOファイルを読み込む
	print STDERR "Loading sequence orientations...";
	read(NSO, my $seq_ori, -s "$input_prefix.nso");
	print STDERR "completed\n";
	
	# NSOファイルを閉じる
	close(NSO);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsi, $input_prefix.nso") unless vec($seq_index, 0, 64) == vec($seq_ori, 0, 64);
	
	# Nanoha Sequence Read (NSR) ファイルを開く
	open(NSR, "<", "$input_prefix.nsr") or &exception::error("failed to open file: $input_prefix.nsr");
	
	# NSRファイルを先頭から8バイト読み込む
	read(NSR, $file_header, 8);
	
	# ファイルヘッダーの一致を確認
	&exception::error("file headers not matched: $input_prefix.nsi, $input_prefix.nsr") unless vec($seq_index, 0, 64) == vec($file_header, 0, 64);
	
	# 変数を宣言
	my @worker_thread = ();
	my $num_error_threads = 0;
	
	# 入出力キューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 指定されたワーカースレッド数で並列処理
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_thread[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューを読み込みながら処理
			while (defined(my $dat = $input->dequeue)) {
				# 変数を宣言
				my %read_seqs = ();
				
				# Nanoha Sequence Read (NSR) ファイルを開く
				open(NSR, "<", "$input_prefix.nsr") or &exception::error("failed to open file: $input_prefix.nsr");
				
				# 各リードについて処理
				foreach my $read_id (unpack("N*", $dat)) {
					# シーケンスリードインデックスを取得
					my $seq_read_index = vec($seq_index, ($read_id - 1) * 2 + 1, 64);
					
					# NSRファイルのファイルポインタをシーケンスリードインデックスの位置にセット
					seek(NSR, $seq_read_index, 0);
					
					# NSRファイルから4バイト読み込む
					read(NSR, my $read_len, 4);
					
					# NSRファイルからリードシーケンスをバイナリ形式で読み込む
					read(NSR, my $read_seq, vec($seq_index, $read_id * 2 + 1, 64) - $seq_read_index - 4);
					
					# リードシーケンスを復号
					$read_seq = join("", map {$convert->to_base(vec($read_seq, $_, 2))} 0..vec($read_len, 0, 32) - 1);
					
					# 復号したシーケンスを追加 (逆鎖の場合は相補鎖に変換)
					$read_seqs{$read_id} = vec($seq_ori, $read_id + 64, 1) ? common::complementary($read_seq) : $read_seq;
				}
				
				# リード長の最小値と最大値を取得
				my ($shortest, $longest) = List::MoreUtils::minmax(map {length($_)} values(%read_seqs));
				
				# 順鎖リード数を算出
				my $num_plus_reads = grep {!vec($seq_ori, $_ + 64, 1)} keys(%read_seqs);
				
				# 逆鎖リード数を算出
				my $num_minus_reads = grep {vec($seq_ori, $_ + 64, 1)} keys(%read_seqs);
				
				# ストランドバイアス値を算出
				my $strand_bias = 1 - ($num_plus_reads < $num_minus_reads ? $num_plus_reads / $num_minus_reads : $num_minus_reads / $num_plus_reads);
				
				# コンセンサスリードシーケンスとして最長のリードシーケンスを取得
				my $consensus_read_seq = [map {$read_seqs{$_}} sort {$a <=> $b} grep {length($read_seqs{$_}) == $longest} keys(%read_seqs)]->[0];
				
				# リード長の長いものから指定した数のリードシーケンスを用いてコンセンサスリードシーケンスを構築 (-c指定時)
				$consensus_read_seq = &main::create_consensus_sequence([map {$read_seqs{$_}} sort {length($read_seqs{$b}) <=> length($read_seqs{$a}) || $a <=> $b} keys(%read_seqs)], $opt{"a"}, $align_method{$opt{"c"}}, $opt{"m"}, -$opt{"n"}, -$opt{"g"}, -$opt{"e"}) if $opt{"c"};
				
				# コンセンサスリードシーケンスからギャップを除去
				$consensus_read_seq =~ s/-//g;
				
				# 出力キューにコンセンサスリード長と最小リード長及び最大リード長、深度、ストランドバイアス値、コンセンサスリードシーケンスをバイナリ形式で追加
				$output->enqueue(pack("NNNNFA*", length($consensus_read_seq), $shortest, $longest, scalar(keys(%read_seqs)), $strand_bias, $consensus_read_seq));
			}
			
			# NSRファイルを閉じる
			close(NSR);
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# データストリームスレッドを作成
	## ここからデータストリームスレッドの処理 ##
	$worker_thread[$opt{"p"}] = threads::async {
		# このスレッドのみを終了可能に変更
		threads->set_thread_exit_only(1);
		
		# 変数を宣言
		my $read_id = 0;
		
		# 出力キューを読み込みながら処理
		while (defined(my $dat = $output->dequeue)) {
			# リードIDを更新
			$read_id++;
			
			# コンセンサスリード長と最小リード長及び最大リード長、深度、ストランドバイアス値、コンセンサスリードシーケンスをデシリアライズして取得
			my ($read_len, $shortest, $longest, $depth, $strand_bias, $read_seq) = unpack("NNNNFA*", $dat);
			
			# シーケンスをFASTA形式で出力
			print ">${input_prefix}_consensus_${read_id} length=$read_len shortest=$shortest longest=$longest depth=$depth strand_bias=$strand_bias\n$read_seq\n";
		}
		
		# スレッドを終了
		return(1);
	};
	## ここまでデータストリームスレッドの処理 ##
	
	# 変数を宣言
	my @clusters = ();
	
	# NSCファイルを読み込みながら処理
	print STDERR "Loading sequence clusters...";
	until (eof(NSC)) {
		# 変数を宣言
		my %cluster = ();
		
		# NSCファイルから4バイト読み込む
		read(NSC, my $num_reads, 4);

		# NSCファイルから所属リード数×4バイト読み込む
		read(NSC, my $read_ids, vec($num_reads, 0, 32) * 4);
		
		# 各所属リードについて処理
		foreach my $read_id (unpack("N*", $read_ids)) {
			# NSRファイルのファイルポインタをシーケンスインデックスの位置にセット
			seek(NSR, vec($seq_index, $read_id * 2 + 1, 64), 0);
			
			# NSRファイルから4バイト読み込む
			read(NSR, my $read_len, 4);

			# クラスターにリードIDとリード長を登録
			$cluster{$read_id} = vec($read_len, 0, 32);
		}
		
		# クラスターリストにクラスターを追加
		push(@clusters, \%cluster);
	}
	print STDERR "completed\n";
	
	# NSRファイルを閉じる
	close(NSR);
	
	# NSCファイルを閉じる
	close(NSC);
	
	# ストランド特異性フラグを取得
	my $strand_speicificity = vec($seq_ori, 64, 1);
	
	# 変数を宣言
	my %binomial_coefficient = ();
	my @p_values = ();
	
	# 各クラスターについて処理
	print STDERR $opt{"c"} ? "Generating a consensus" : "Extracting the longest", " sequence read from each cluster...";
	foreach my $cluster (@clusters) {
		# ストランド特異性フラグが立っている場合はp-valueに0を代入して以下の処理をスキップ
		push(@p_values, 0) and next if $strand_speicificity;
		
		# 変数を宣言
		my @num_reads = (0, 0);
		my $num_reads = keys(%{$cluster});
		
		# シーケンス方向ごとにリード数を集計
		$num_reads[vec($seq_ori, $_ + 64, 1)]++ foreach keys(%{$cluster});
		
		# 二項係数を算出
		&calc_binomial_coefficient(\%binomial_coefficient, $num_reads, List::Util::min(@num_reads));
		
		# ストランドバイアスのp-valueを算出
		push(@p_values, List::Util::sum(map {$binomial_coefficient{$num_reads}->[$_]} 0..List::Util::min(@num_reads)) / 2 ** $num_reads);
	}
	
	# 変数を宣言
	my @q_values = ();
	my @processed_clusters = ();
	
	# 各クラスターについてストランドバイアスのq-valueを算出
	push(@processed_clusters, $_) and $q_values[$_] = $p_values[$_] * @clusters / @processed_clusters foreach sort {$p_values[$a] <=> $p_values[$b] || $a <=> $b} 0..$#clusters;
	@q_values[reverse(@processed_clusters)] = map {$q_values[$_]} List::Util::reductions {$q_values[$a] < $q_values[$b] ? $a : $b} reverse(@processed_clusters);
	
	# ストランドバイアスのq-valueが指定値以上かつクラスターサイズが指定値以上の対象クラスターを取得
	my @target_clusters = grep {keys(%{$_}) >= $opt{"d"}} map {$clusters[$_]} grep {$strand_speicificity or $q_values[$_] >= $opt{"q"}} 0..$#clusters;
	
	# 各対象クラスターについて実行中のスレッド数が指定値+1と一致している場合は所属リードのリードIDをシリアライズして入力キューに追加
	threads->list(threads::running) == $opt{"p"} + 1 and $input->enqueue(pack("N*", keys(%{$_}))) foreach @target_clusters;
	
	# 入力キューを終了
	$input->end;
	
	# 並列処理の各ワーカースレッドが終了するまで待機
	$worker_thread[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 出力キューを終了
	$output->end;
	
	# データストリームスレッドが終了するまで待機
	$worker_thread[$opt{"p"}]->join or &exception::error("data stream thread abnormally exited");
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	return(1);
}

# 二項係数を算出 correct::calc_binomial_coefficient(二項係数ハッシュリファレンス, 整数, 整数)
sub calc_binomial_coefficient {
	# 引数を取得
	my ($binomial_coefficient, $n, $r) = @_;
	
	# 初期値を定義
	$binomial_coefficient->{$n}->[0] = Math::BigFloat->new(1) unless exists($binomial_coefficient->{$n});
	
	# 未算出の部分を算出
	$binomial_coefficient->{$n}->[$_] = $binomial_coefficient->{$n}->[$_ - 1] * ($n - $_ + 1) / $_ foreach scalar(@{$binomial_coefficient->{$n}})..$r;
	return(1);
}

## ここから共通処理のパッケージ ##
package common;

# ファイル確認 common::check_files(ファイルリストリファレンス)
sub check_files {
	# 引数を取得
	my ($file_list) = @_;
	
	# 各ファイルについて処理
	foreach my $file (@{$file_list}) {
		next if !$file;
		&exception::error("file not found: $file") if !-f $file;
		&exception::error("file unreadable: $file") if !-r $file;
		&exception::error("null file specified: $file") if !-s $file;
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

# サブルーチンを追加

# パッケージを追加
### 編集範囲 終了 ###
__END__
__C__
// LLCSに基づくJaccard類似度を算出 classify::calc_LLCS_Jaccard_similarities_0(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性フラグ)
void calc_LLCS_Jaccard_similarities_0(AV* read_seqs, AV* read_lens, unsigned int strand_speicificity) {
	// 戻り値としてスタックの利用を宣言
	Inline_Stack_Vars;
	
	// スタックをリセット
	Inline_Stack_Reset;
	
	// クエリーリードシーケンスを取得
	const uint64_t* query_seq = SvPV_nolen(*av_fetch(read_seqs, 0, 0));
	
	// クエリーリード長を取得
	const uint32_t query_len = SvUV(*av_fetch(read_lens, 0, 0));
	
	// クエリーブロック数を算出
	const uint32_t num_query_blocks = (query_len >> 6) + ((query_len & 0x3F) > 0);
	
	// 変数を宣言
	uint64_t* query_matrix[4];
	uint64_t query_base[2];
	uint64_t mask = 0xFFFFFFFFFFFFFFFF;
	uint64_t for_llcs;
	uint64_t rev_llcs;
	uint64_t* v;
	uint64_t u;
	uint64_t w;
	uint64_t x;
	
	// メモリを確保
	query_matrix[0] = malloc(sizeof(uint64_t) * num_query_blocks);
	query_matrix[1] = malloc(sizeof(uint64_t) * num_query_blocks);
	query_matrix[2] = malloc(sizeof(uint64_t) * num_query_blocks);
	query_matrix[3] = malloc(sizeof(uint64_t) * num_query_blocks);
	v = malloc(sizeof(uint64_t) * (num_query_blocks + 3));
	
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
		
		// 計算要素を初期化
		for_llcs = 0;
		for (uint32_t i = 0;i < num_query_blocks + 3;i++) {v[i] = 0xFFFFFFFFFFFFFFFF;}
		
		// 順鎖に対するLLCSを算出
		for (uint32_t i = 0;i < num_blocks;i++) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = 0;j < end;j += 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint64_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = v[k] & query_matrix[base][k];
					w = v[k] + u + carry;
					carry = w < v[k];
					v[k] = w | (v[k] - u);
				}
			}
		}
		for (uint32_t i = 0;i < num_query_blocks;i += 4) {
			x = 0;
			for (uint32_t j = 0;j < 4;j += 2) {
				w = 0;
				for (uint32_t k = 0;k < 2;k++) {
					u = ~v[i + j + k];
					u = (u & 0x5555555555555555) + (u >> 1 & 0x5555555555555555);
					u = (u & 0x3333333333333333) + (u >> 2 & 0x3333333333333333);
					w += u;
				}
				w = (w & 0x0F0F0F0F0F0F0F0F) + (w >> 4 & 0x0F0F0F0F0F0F0F0F);
				w = (w & 0x00FF00FF00FF00FF) + (w >> 8 & 0x00FF00FF00FF00FF);
				x += w;
			}
			x = (x & 0x0000FFFF0000FFFF) + (x >> 16 & 0x0000FFFF0000FFFF);
			x = (x & 0x00000000FFFFFFFF) + (x >> 32 & 0x00000000FFFFFFFF);
			for_llcs += x;
		}
		
		// ストランド特異性フラグが立っている場合は順鎖に対するLLCSに基づくJaccard類似度を算出して以下の処理をスキップ
		if (strand_speicificity) {
			Inline_Stack_Push(newSVnv((double)for_llcs / (query_len + read_len - for_llcs)));
			continue;
		}
		
		// リードシーケンスを相補的な塩基に変換
		for (uint32_t i = 0;i < num_blocks;i++) {read_seq[i] = ~read_seq[i];}
		
		// 計算要素を初期化
		rev_llcs = 0;
		for (uint32_t i = 0;i < num_query_blocks + 3;i++) {v[i] = 0xFFFFFFFFFFFFFFFF;}
		
		// 逆鎖に対するLLCSを算出
		for (uint32_t i = num_blocks - 1;i < num_blocks;i--) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = end - 2;j < end;j -= 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint32_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = v[k] & query_matrix[base][k];
					w = v[k] + u + carry;
					carry = w < v[k];
					v[k] = w | (v[k] - u);
				}
			}
		}
		for (uint32_t i = 0;i < num_query_blocks;i += 4) {
			x = 0;
			for (uint32_t j = 0;j < 4;j += 2) {
				w = 0;
				for (uint32_t k = 0;k < 2;k++) {
					u = ~v[i + j + k];
					u = (u & 0x5555555555555555) + (u >> 1 & 0x5555555555555555);
					u = (u & 0x3333333333333333) + (u >> 2 & 0x3333333333333333);
					w += u;
				}
				w = (w & 0x0F0F0F0F0F0F0F0F) + (w >> 4 & 0x0F0F0F0F0F0F0F0F);
				w = (w & 0x00FF00FF00FF00FF) + (w >> 8 & 0x00FF00FF00FF00FF);
				x += w;
			}
			x = (x & 0x0000FFFF0000FFFF) + (x >> 16 & 0x0000FFFF0000FFFF);
			x = (x & 0x00000000FFFFFFFF) + (x >> 32 & 0x00000000FFFFFFFF);
			rev_llcs += x;
		}
		
		// LLCSに基づくJaccard類似度を算出 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		Inline_Stack_Push(newSVnv(for_llcs >= rev_llcs ? (double)for_llcs / (query_len + read_len - for_llcs) : -(double)rev_llcs / (query_len + read_len - rev_llcs)));
	}
	
	CLEANUP:
	// メモリを解放
	free(query_matrix[0]);
	free(query_matrix[1]);
	free(query_matrix[2]);
	free(query_matrix[3]);
	free(v);
	
	// スタックを終了
	Inline_Stack_Done;
	
	// LLCSに基づくJaccard類似度を返す
	return;
}

__C__
#define word_size_index 1
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

// LLCSに基づくJaccard類似度を算出 (SSE4.1) main::calc_LLCS_Jaccard_similarities_1(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性フラグ)
void calc_LLCS_Jaccard_similarities_1(AV* read_seqs, AV* read_lens, unsigned int strand_speicificity) {
	// 戻り値としてスタックの利用を宣言
	Inline_Stack_Vars;
	
	// スタックをリセット
	Inline_Stack_Reset;
	
	// クエリーリードシーケンスを取得
	const uint64_t* query_seq = SvPV_nolen(*av_fetch(read_seqs, 0, 0));
	
	// クエリーリード長を取得
	const uint32_t query_len = SvUV(*av_fetch(read_lens, 0, 0));
	
	// クエリーブロック数を算出
	const uint32_t num_query_blocks = (query_len >> word_size_index + 6) + ((query_len & ((0x40 << word_size_index) - 1)) > 0);
	
	// 定数を定義
	const vector_int mask_ff = vec_splat_qword(0xFFFFFFFFFFFFFFFF);
	const vector_int mask_55 = vec_splat_qword(0x5555555555555555);
	const vector_int mask_0f = vec_splat_qword(0x0F0F0F0F0F0F0F0F);
	const vector_int compress_table1 = vec_set_qword(0x0F0E0B0A0D0C0908, 0x0706030205040100);
	const vector_int compress_table2 = vec_set_qword(0x1001100110011001, 0x1001100110011001);
	const vector_int maskend_table1 = vec_set_qword(0xF1F2F3F4F5F6F7F8, 0xF9FAFBFCFDFEFF00);
	const vector_int maskend_table2 = vec_set_qword(0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF);
	const vector_int extract_table = vec_set_qword(0x0000000000000002, 0x0000000000000001);
	const vector_int popcnt_table1 = vec_set_qword(0x0403030203020201, 0x0302020102010100);
	const vector_int popcnt_table2 = vec_set_qword(0x0405050605060607, 0x0506060706070708);
	
	// 変数を宣言
	vector_int* query_matrix[4];
	vector_int query_base[2];
	vector_int mask = mask_ff;
	vector_int for_llcs;
	vector_int rev_llcs;
	vector_int* v;
	vector_int u;
	vector_int w;
	vector_int x;
	
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
		
		// 計算要素を初期化
		for_llcs = vec_zero;
		for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = mask_ff;}
		
		// 順鎖に対するLLCSを算出
		for (uint32_t i = 0;i < num_blocks;i++) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = 0;j < end;j += 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint64_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = vec_and(v[k], query_matrix[base][k]);
					w = vec_add_qword(v[k], u);
					carry = vec_movemask_double((vector_double)vec_cmpeq_qword(w, mask_ff)) + (vec_movemask_double((vector_double)vec_or(vec_andnot(w, v[k]), u)) << 1) + (carry >> (0x01 << word_size_index));
					v[k] = vec_or(vec_sub_qword(w, vec_cmpeq_qword(vec_and(vec_splat_qword(carry), extract_table), extract_table)), vec_xor(v[k], u));
				}
			}
		}
		u = num_query_blocks & 0x01 ? v[0] : mask_ff;
		for (uint32_t i = num_query_blocks & 0x01;i < num_query_blocks;i += 2) {
			x = vec_xor(v[i], v[i + 1]);
			w = vec_or(vec_and(v[i], v[i + 1]), vec_and(u, x));
			u = vec_xor(u, x);
			for_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(w, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(w, 4), mask_0f))), for_llcs);
		}
		for_llcs = vec_add_qword(for_llcs, for_llcs);
		for_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(u, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(u, 4), mask_0f))), for_llcs);
		
		// ストランド特異性フラグが立っている場合は順鎖に対するLLCSに基づくJaccard類似度を算出して以下の処理をスキップ
		if (strand_speicificity) {
			for_llcs = vec_add_dword(vec_shuffle_dword(for_llcs, 0x4E), for_llcs);
			uint32_t llcs = vec_cast_dword(for_llcs);
			Inline_Stack_Push(newSVnv((double)llcs / (query_len + read_len - llcs)));
			continue;
		}
		
		// リードシーケンスを相補的な塩基に変換
		for (uint32_t i = 0;i < num_blocks;i++) {read_seq[i] = ~read_seq[i];}
		
		// 計算要素を初期化
		rev_llcs = vec_zero;
		for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = mask_ff;}
		
		// 逆鎖に対するLLCSを算出
		for (uint32_t i = num_blocks - 1;i < num_blocks;i--) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = end - 2;j < end;j -= 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint32_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = vec_and(v[k], query_matrix[base][k]);
					w = vec_add_qword(v[k], u);
					carry = vec_movemask_double((vector_double)vec_cmpeq_qword(w, mask_ff)) + (vec_movemask_double((vector_double)vec_or(vec_andnot(w, v[k]), u)) << 1) + (carry >> (0x01 << word_size_index));
					v[k] = vec_or(vec_sub_qword(w, vec_cmpeq_qword(vec_and(vec_splat_qword(carry), extract_table), extract_table)), vec_xor(v[k], u));
				}
			}
		}
		u = num_query_blocks & 0x01 ? v[0] : mask_ff;
		for (uint32_t i = num_query_blocks & 0x01;i < num_query_blocks;i += 2) {
			x = vec_xor(v[i], v[i + 1]);
			w = vec_or(vec_and(v[i], v[i + 1]), vec_and(u, x));
			u = vec_xor(u, x);
			rev_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(w, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(w, 4), mask_0f))), rev_llcs);
		}
		rev_llcs = vec_add_qword(rev_llcs, rev_llcs);
		rev_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(u, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(u, 4), mask_0f))), rev_llcs);
		
		// LLCSに基づくJaccard類似度を算出 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		w = vec_add_dword(for_llcs, rev_llcs);
		w = vec_add_dword(vec_shuffle_dword(w, 0x4E), w);
		u = vec_sub_dword(for_llcs, rev_llcs);
		u = vec_add_dword(vec_shuffle_dword(u, 0x4E), u);
		x = vec_cmpgt_dword(vec_zero, u);
		int32_t llcs = vec_cast_dword(vec_srai_dword(vec_add_dword(vec_sub_dword(vec_xor(w, x), x), u), 1));
		Inline_Stack_Push(newSVnv((double)llcs / (query_len + read_len - abs(llcs))));
	}
	
	CLEANUP:
	// メモリを解放
	_mm_free(query_matrix[0]);
	_mm_free(query_matrix[1]);
	_mm_free(query_matrix[2]);
	_mm_free(query_matrix[3]);
	_mm_free(v);
	
	// スタックを終了
	Inline_Stack_Done;
	
	// LLCSに基づくJaccard類似度を返す
	return;
}

__C__
#define word_size_index 2
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

// LLCSに基づくJaccard類似度を算出 (AVX2) main::calc_LLCS_Jaccard_similarities_2(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性フラグ)
void calc_LLCS_Jaccard_similarities_2(AV* read_seqs, AV* read_lens, unsigned int strand_speicificity) {
	// 戻り値としてスタックの利用を宣言
	Inline_Stack_Vars;
	
	// スタックをリセット
	Inline_Stack_Reset;
	
	// クエリーリードシーケンスを取得
	const uint64_t* query_seq = SvPV_nolen(*av_fetch(read_seqs, 0, 0));
	
	// クエリーリード長を取得
	const uint32_t query_len = SvUV(*av_fetch(read_lens, 0, 0));
	
	// クエリーブロック数を算出
	const uint32_t num_query_blocks = (query_len >> word_size_index + 6) + ((query_len & ((0x40 << word_size_index) - 1)) > 0);
	
	// 定数を定義
	const vector_int mask_ff = vec_splat_qword(0xFFFFFFFFFFFFFFFF);
	const vector_int mask_55 = vec_splat_qword(0x5555555555555555);
	const vector_int mask_0f = vec_splat_qword(0x0F0F0F0F0F0F0F0F);
	const vector_int compress_table1 = vec_set_qword(0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100);
	const vector_int compress_table2 = vec_set_qword(0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001);
	const vector_int maskend_table1 = vec_set_qword(0xE1E2E3E4E5E6E7E8, 0xE9EAEBECEDEEEFF0, 0xF1F2F3F4F5F6F7F8, 0xF9FAFBFCFDFEFF00);
	const vector_int maskend_table2 = vec_set_qword(0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF);
	const vector_int extract_table = vec_set_qword(0x0000000000000008, 0x0000000000000004, 0x0000000000000002, 0x0000000000000001);
	const vector_int popcnt_table1 = vec_set_qword(0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100);
	const vector_int popcnt_table2 = vec_set_qword(0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708);
	const vector_int permute_table = vec_set_qword(0x0000000300000002, 0x0000000100000000, 0x0000000700000006, 0x0000000500000004);
	
	// 変数を宣言
	vector_int* query_matrix[4];
	vector_int query_base[2];
	vector_int mask = mask_ff;
	vector_int for_llcs;
	vector_int rev_llcs;
	vector_int* v;
	vector_int u;
	vector_int w;
	vector_int x;
	
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
		
		// 計算要素を初期化
		for_llcs = vec_zero;
		for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = mask_ff;}
		
		// 順鎖に対するLLCSを算出
		for (uint32_t i = 0;i < num_blocks;i++) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = 0;j < end;j += 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint64_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = vec_and(v[k], query_matrix[base][k]);
					w = vec_add_qword(v[k], u);
					carry = vec_movemask_double((vector_double)vec_cmpeq_qword(w, mask_ff)) + (vec_movemask_double((vector_double)vec_or(vec_andnot(w, v[k]), u)) << 1) + (carry >> (0x01 << word_size_index));
					v[k] = vec_or(vec_sub_qword(w, vec_cmpeq_qword(vec_and(vec_splat_qword(carry), extract_table), extract_table)), vec_xor(v[k], u));
				}
			}
		}
		u = num_query_blocks & 0x01 ? v[0] : mask_ff;
		for (uint32_t i = num_query_blocks & 0x01;i < num_query_blocks;i += 2) {
			x = vec_xor(v[i], v[i + 1]);
			w = vec_or(vec_and(v[i], v[i + 1]), vec_and(u, x));
			u = vec_xor(u, x);
			for_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(w, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(w, 4), mask_0f))), for_llcs);
		}
		for_llcs = vec_add_qword(for_llcs, for_llcs);
		for_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(u, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(u, 4), mask_0f))), for_llcs);
		
		// ストランド特異性フラグが立っている場合は順鎖に対するLLCSに基づくJaccard類似度を算出して以下の処理をスキップ
		if (strand_speicificity) {
			for_llcs = vec_add_dword(vec_shuffle_dword(for_llcs, 0x4E), for_llcs);
			for_llcs = vec_add_dword(vec_permute_dword(for_llcs, permute_table), for_llcs);
			uint32_t llcs = vec_cast_dword(for_llcs);
			Inline_Stack_Push(newSVnv((double)llcs / (query_len + read_len - llcs)));
			continue;
		}
		
		// リードシーケンスを相補的な塩基に変換
		for (uint32_t i = 0;i < num_blocks;i++) {read_seq[i] = ~read_seq[i];}
		
		// 計算要素を初期化
		rev_llcs = vec_zero;
		for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = mask_ff;}
		
		// 逆鎖に対するLLCSを算出
		for (uint32_t i = num_blocks - 1;i < num_blocks;i--) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = end - 2;j < end;j -= 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint32_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = vec_and(v[k], query_matrix[base][k]);
					w = vec_add_qword(v[k], u);
					carry = vec_movemask_double((vector_double)vec_cmpeq_qword(w, mask_ff)) + (vec_movemask_double((vector_double)vec_or(vec_andnot(w, v[k]), u)) << 1) + (carry >> (0x01 << word_size_index));
					v[k] = vec_or(vec_sub_qword(w, vec_cmpeq_qword(vec_and(vec_splat_qword(carry), extract_table), extract_table)), vec_xor(v[k], u));
				}
			}
		}
		u = num_query_blocks & 0x01 ? v[0] : mask_ff;
		for (uint32_t i = num_query_blocks & 0x01;i < num_query_blocks;i += 2) {
			x = vec_xor(v[i], v[i + 1]);
			w = vec_or(vec_and(v[i], v[i + 1]), vec_and(u, x));
			u = vec_xor(u, x);
			rev_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(w, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(w, 4), mask_0f))), rev_llcs);
		}
		rev_llcs = vec_add_qword(rev_llcs, rev_llcs);
		rev_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(u, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(u, 4), mask_0f))), rev_llcs);
		
		// LLCSに基づくJaccard類似度を算出 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		w = vec_add_dword(for_llcs, rev_llcs);
		w = vec_add_dword(vec_shuffle_dword(w, 0x4E), w);
		w = vec_add_dword(vec_permute_dword(w, permute_table), w);
		u = vec_sub_dword(for_llcs, rev_llcs);
		u = vec_add_dword(vec_shuffle_dword(u, 0x4E), u);
		u = vec_add_dword(vec_permute_dword(u, permute_table), u);
		x = vec_cmpgt_dword(vec_zero, u);
		int32_t llcs = vec_cast_dword(vec_srai_dword(vec_add_dword(vec_sub_dword(vec_xor(w, x), x), u), 1));
		Inline_Stack_Push(newSVnv((double)llcs / (query_len + read_len - abs(llcs))));
	}
	
	CLEANUP:
	// メモリを解放
	_mm_free(query_matrix[0]);
	_mm_free(query_matrix[1]);
	_mm_free(query_matrix[2]);
	_mm_free(query_matrix[3]);
	_mm_free(v);
	
	// スタックを終了
	Inline_Stack_Done;
	
	// LLCSに基づくJaccard類似度を返す
	return;
}

__C__
#define word_size_index 3
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

// LLCSに基づくJaccard類似度を算出 (AVX-512) main::calc_LLCS_Jaccard_similarities_3(リードシーケンスリストリファレンス, リード長リストリファレンス, ストランド特異性フラグ)
void calc_LLCS_Jaccard_similarities_3(AV* read_seqs, AV* read_lens, unsigned int strand_speicificity) {
	// 戻り値としてスタックの利用を宣言
	Inline_Stack_Vars;
	
	// スタックをリセット
	Inline_Stack_Reset;
	
	// クエリーリードシーケンスを取得
	const uint64_t* query_seq = SvPV_nolen(*av_fetch(read_seqs, 0, 0));
	
	// クエリーリード長を取得
	const uint32_t query_len = SvUV(*av_fetch(read_lens, 0, 0));
	
	// クエリーブロック数を算出
	const uint32_t num_query_blocks = (query_len >> word_size_index + 6) + ((query_len & ((0x40 << word_size_index) - 1)) > 0);
	
	// 定数を定義
	const vector_int mask_ff = vec_splat_qword(0xFFFFFFFFFFFFFFFF);
	const vector_int mask_55 = vec_splat_qword(0x5555555555555555);
	const vector_int mask_0f = vec_splat_qword(0x0F0F0F0F0F0F0F0F);
	const vector_int compress_table1 = vec_set_qword(0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100, 0x0F0E0B0A0D0C0908, 0x0706030205040100);
	const vector_int compress_table2 = vec_set_qword(0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001, 0x1001100110011001);
	const vector_int maskend_table1 = vec_set_qword(0xC1C2C3C4C5C6C7C8, 0xC9CACBCCCDCECFD0, 0xD1D2D3D4D5D6D7D8, 0xD9DADBDCDDDEDFE0, 0xE1E2E3E4E5E6E7E8, 0xE9EAEBECEDEEEFF0, 0xF1F2F3F4F5F6F7F8, 0xF9FAFBFCFDFEFF00);
	const vector_int maskend_table2 = vec_set_qword(0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF, 0x7F3F1F0F070301FF, 0xFFFFFFFFFFFFFFFF);
	const vector_int popcnt_table1 = vec_set_qword(0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100, 0x0403030203020201, 0x0302020102010100);
	const vector_int popcnt_table2 = vec_set_qword(0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708, 0x0405050605060607, 0x0506060706070708);
	const vector_int permute_table = vec_set_qword(0x0000000700000006, 0x0000000500000004, 0x0000000300000002, 0x0000000100000000, 0x0000000F0000000E, 0x0000000D0000000C, 0x0000000B0000000A, 0x0000000900000008);
	
	// 変数を宣言
	vector_int* query_matrix[4];
	vector_int query_base[2];
	vector_int mask = mask_ff;
	vector_int for_llcs;
	vector_int rev_llcs;
	vector_int* v;
	vector_int u;
	vector_int w;
	
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
		
		// 計算要素を初期化
		for_llcs = vec_zero;
		for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = mask_ff;}
		
		// 順鎖に対するLLCSを算出
		for (uint32_t i = 0;i < num_blocks;i++) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = 0;j < end;j += 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint64_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = vec_and(v[k], query_matrix[base][k]);
					w = vec_add_qword(v[k], u);
					carry = vec_cmpeq_qword(w, mask_ff) + (vec_cmplt_qword(w, v[k]) << 1) + (carry >> (0x01 << word_size_index));
					v[k] = vec_or(vec_mask_sub_qword(w, carry, w, mask_ff), vec_xor(v[k], u));
				}
			}
		}
		u = num_query_blocks & 0x01 ? v[0] : mask_ff;
		for (uint32_t i = num_query_blocks & 0x01;i < num_query_blocks;i += 2) {
			w = vec_ternarylogic(u, v[i], v[i + 1], 0xE8);
			u = vec_ternarylogic(u, v[i], v[i + 1], 0x96);
			for_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(w, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(w, 4), mask_0f))), for_llcs);
		}
		for_llcs = vec_add_qword(for_llcs, for_llcs);
		for_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(u, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(u, 4), mask_0f))), for_llcs);
		
		// ストランド特異性フラグが立っている場合は順鎖に対するLLCSに基づくJaccard類似度を算出して以下の処理をスキップ
		if (strand_speicificity) {
			for_llcs = vec_add_dword(vec_shuffle_dword(for_llcs, 0x4E), for_llcs);
			for_llcs = vec_add_dword(vec_permute_qword(for_llcs, 0x4E), for_llcs);
			for_llcs = vec_add_dword(vec_permute_dword(for_llcs, permute_table), for_llcs);
			uint32_t llcs = vec_cast_dword(for_llcs);
			Inline_Stack_Push(newSVnv((double)llcs / (query_len + read_len - llcs)));
			continue;
		}
		
		// リードシーケンスを相補的な塩基に変換
		for (uint32_t i = 0;i < num_blocks;i++) {read_seq[i] = ~read_seq[i];}
		
		// 計算要素を初期化
		rev_llcs = vec_zero;
		for (uint32_t i = 0;i < num_query_blocks;i++) {v[i] = mask_ff;}
		
		// 逆鎖に対するLLCSを算出
		for (uint32_t i = num_blocks - 1;i < num_blocks;i--) {
			uint32_t end = i + !!(read_len & 0x1F) < num_blocks ? 64 : read_len << 1 & 0x3E;
			for (uint32_t j = end - 2;j < end;j -= 2) {
				uint32_t base = read_seq[i] >> j & 0x03;
				uint32_t carry = 0;
				for (uint32_t k = 0;k < num_query_blocks;k++) {
					u = vec_and(v[k], query_matrix[base][k]);
					w = vec_add_qword(v[k], u);
					carry = vec_cmpeq_qword(w, mask_ff) + (vec_cmplt_qword(w, v[k]) << 1) + (carry >> (0x01 << word_size_index));
					v[k] = vec_or(vec_mask_sub_qword(w, carry, w, mask_ff), vec_xor(v[k], u));
				}
			}
		}
		u = num_query_blocks & 0x01 ? v[0] : mask_ff;
		for (uint32_t i = num_query_blocks & 0x01;i < num_query_blocks;i += 2) {
			w = vec_ternarylogic(u, v[i], v[i + 1], 0xE8);
			u = vec_ternarylogic(u, v[i], v[i + 1], 0x96);
			rev_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(w, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(w, 4), mask_0f))), rev_llcs);
		}
		rev_llcs = vec_add_qword(rev_llcs, rev_llcs);
		rev_llcs = vec_add_qword(vec_sad(vec_shuffle_byte(popcnt_table1, vec_and(u, mask_0f)), vec_shuffle_byte(popcnt_table2, vec_and(vec_srli_qword(u, 4), mask_0f))), rev_llcs);
		
		// LLCSに基づくJaccard類似度を算出 (逆鎖に対するLLCSが順鎖に対するLLCSより大きい場合は負値にする)
		w = vec_add_dword(for_llcs, rev_llcs);
		w = vec_add_dword(vec_shuffle_dword(w, 0x4E), w);
		w = vec_add_dword(vec_permute_qword(w, 0x4E), w);
		w = vec_add_dword(vec_permute_dword(w, permute_table), w);
		u = vec_sub_dword(for_llcs, rev_llcs);
		u = vec_add_dword(vec_shuffle_dword(u, 0x4E), u);
		u = vec_add_dword(vec_permute_qword(u, 0x4E), u);
		u = vec_add_dword(vec_permute_dword(u, permute_table), u);
		int32_t llcs = vec_cast_dword(vec_srai_dword(vec_add_dword(vec_mask_sub_dword(w, vec_cmpgt_dword(vec_zero, u), vec_xor(w, mask_ff), mask_ff), u), 1));
		Inline_Stack_Push(newSVnv((double)llcs / (query_len + read_len - abs(llcs))));
	}
	
	CLEANUP:
	// メモリを解放
	_mm_free(query_matrix[0]);
	_mm_free(query_matrix[1]);
	_mm_free(query_matrix[2]);
	_mm_free(query_matrix[3]);
	_mm_free(v);
	
	// スタックを終了
	Inline_Stack_Done;
	
	// LLCSに基づくJaccard類似度を返す
	return;
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
	for (const auto& aligned_seq: aligned_read_seqs) {
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
