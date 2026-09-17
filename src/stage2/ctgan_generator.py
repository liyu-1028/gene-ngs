# --- src/stage2/ctgan_generator.py (v2, MSK-IMPACT 50K 量级) ---
import pandas as pd
from sdv.single_table import CTGANSynthesizer
from sdv.metadata import SingleTableMetadata

class ClinicalTabularGenerator:
    """
    CTGAN 临床表格生成器。
    v2: 超参数按训练集规模自适应——
      - 小样本 (n<100, 旧 24 例时代): epochs=300, batch_size=10
      - 大样本 (n>1e4, MSK-IMPACT 50K 时代): epochs=100, batch_size≈N/30
        总更新步数 ≈ 3,000 量级，CPU 训练 ~30-60 分钟
    """

    def __init__(self, epochs: int = None, batch_size: int = None):
        self.epochs = epochs
        self.batch_size = batch_size
        self.model = None
        self.metadata = None

    def fit(self, real_data_df: pd.DataFrame):
        """基于真实数据训练 CTGAN"""
        n = len(real_data_df)
        epochs = self.epochs or (300 if n < 100 else 100)
        if self.batch_size:
            batch_size = self.batch_size
        else:
            batch_size = 10 if n < 100 else max(500, min(2000, round(n / 30 / 10) * 10))
        print(f"Training CTGAN (n={n:,}, epochs={epochs}, batch_size={batch_size}) ...")

        # 1. 自动推断数据元数据（列类型、主键等）
        self.metadata = SingleTableMetadata()
        self.metadata.detect_from_dataframe(real_data_df)

        # 强制纠正 SDV 对 PII 的误判，将其设为分类变量
        for col in real_data_df.columns:
            if 'pathological_diagnosis' in col or 'sample_type' in col:
                self.metadata.update_column(column_name=col, sdtype='categorical')

        # 2. 初始化 CTGAN 综合器
        self.model = CTGANSynthesizer(
            self.metadata,
            epochs=epochs,
            batch_size=batch_size,
            generator_dim=(128, 128),
            discriminator_dim=(128, 128),
            verbose=True
        )
        self.model.fit(real_data_df)
        print("Training Complete.")

    def sample(self, num_rows: int) -> pd.DataFrame:
        """生成合成数据"""
        if self.model is None:
            raise ValueError("Model not trained yet.")
        return self.model.sample(num_rows=num_rows)
