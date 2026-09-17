#!/bin/bash
PROJECT_ROOT="${PROJECT_ROOT:-$HOME/gene-ngs}"  # override with env var
# tstr_watch.sh — 监视 50K TSTR (PID 4912) 完成, 自动验证并更新报告
TARGET=${PROJECT_ROOT}/data/synthetic/evaluation/utility_tstr_metrics.json
REPORT=${PROJECT_ROOT}/OVERNIGHT_REPORT.md
while true; do
    if [ -f "$TARGET" ]; then
        echo "═══ 50K TSTR 完成 $(date '+%H:%M') ═══" >> "$REPORT"
        cat "$TARGET" >> "$REPORT"
        echo "" >> "$REPORT"
        # 验证 JSON 可解析
        ${PROJECT_ROOT}/.venv/bin/python -c "
import json
d = json.load(open('$TARGET'))
p = d['Performance']
print('[watcher] 50K TSTR 结果: TRTR AUROC=%s, TSTR AUROC=%s, Δ=%s, SHAP rho=%s' % (
    p['TRTR']['AUROC'], p['TSTR']['AUROC'], p['Delta_AUROC'],
    d['Explainability']['SHAP_Rank_Correlation']))" >> "$REPORT" 2>&1
        exit 0
    fi
    # 进程死亡且无输出 = 异常
    if ! kill -0 4912 2>/dev/null; then
        echo "═══ ⚠ TSTR 进程 4912 已退出但无输出文件 $(date '+%H:%M') ═══" >> "$REPORT"
        exit 1
    fi
    sleep 120
done
