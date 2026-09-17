#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compress_watch_v2.py — 崩溃恢复期连续压缩守护
持续扫描 reads 目录, 将"完成且稳定"的患者 FASTQ 目录压成单个 <PT>.fastq.gz 后删除原文件,
为后续仿真腾出磁盘。安全机制:
  - 目录 mtime 需静默 >= min-age 分钟 (避免与 synggen 写入竞态)
  - 任一 fastq < 1MB 视为不完整, 跳过
  - gzip -t 校验通过才删除原文件
用法: python compress_watch_v2.py --reads-root data/synthetic/reads --reads-root data/synthetic_2017/reads --stop-file /tmp/compress.stop
"""
import argparse, glob, gzip, os, shutil, subprocess, time

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reads-root', action='append', required=True)
    ap.add_argument('--min-age-min', type=float, default=10.0)
    ap.add_argument('--stop-file', default='/tmp/compress.stop')
    ap.add_argument('--poll-sec', type=int, default=30)
    args = ap.parse_args()
    os.nice(10)

    min_age = args.min_age_min * 60
    done = 0
    print(f"[compress-watch] 启动: roots={args.reads_root}, 停止文件={args.stop_file}", flush=True)
    while not os.path.exists(args.stop_file):
        worked = False
        for root in args.reads_root:
            if not os.path.isdir(root):
                continue
            for d in sorted(glob.glob(os.path.join(root, 'virtual_pt_*'))):
                if not os.path.isdir(d):
                    continue
                gz_path = d.rstrip('/') + '.fastq.gz'
                if os.path.exists(gz_path):
                    continue
                fastqs = sorted(glob.glob(os.path.join(d, '*.fastq')))
                if not fastqs:
                    continue
                # 完整性: 所有 fastq >= 1MB, 且目录已静默 min_age
                if any(os.path.getsize(f) < 1_000_000 for f in fastqs):
                    continue
                if time.time() - os.path.getmtime(d) < min_age:
                    continue
                # 压缩 (cat + gzip 流式)
                print(f"[compress-watch] {d} ({len(fastqs)} files) -> {gz_path}", flush=True)
                try:
                    with gzip.open(gz_path + '.tmp', 'wb', compresslevel=6) as out:
                        for fq in fastqs:
                            with open(fq, 'rb') as fh:
                                shutil.copyfileobj(fh, out, length=8 * 1024 * 1024)
                    # 校验
                    r = subprocess.run(['gzip', '-t', gz_path + '.tmp'],
                                       capture_output=True, timeout=600)
                    if r.returncode != 0:
                        print(f"[compress-watch] ❌ gzip -t 失败, 保留原文件: {d}", flush=True)
                        os.remove(gz_path + '.tmp')
                        continue
                    os.rename(gz_path + '.tmp', gz_path)
                    # 删除原文件 (fastq + 可再生成的 extended.bed)
                    for f in glob.glob(os.path.join(d, '*.fastq')) + \
                             glob.glob(os.path.join(d, '*.extended.bed')):
                        os.remove(f)
                    try:
                        os.rmdir(d)
                    except OSError:
                        pass
                    done += 1
                    worked = True
                    print(f"[compress-watch] ✅ 累计 {done} 个患者", flush=True)
                except Exception as e:
                    print(f"[compress-watch] ❌ {d}: {e}", flush=True)
                    if os.path.exists(gz_path + '.tmp'):
                        os.remove(gz_path + '.tmp')
        if not worked:
            time.sleep(args.poll_sec)
    print(f"[compress-watch] 收到停止信号, 共压缩 {done} 个患者, 退出", flush=True)

if __name__ == '__main__':
    main()
