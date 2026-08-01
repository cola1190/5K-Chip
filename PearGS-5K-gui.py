import sys
import os
import glob
import pickle
import numpy as np
import pandas as pd

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QTextEdit, QFileDialog,
    QMessageBox, QFrame, QListWidget, QListWidgetItem
)
from PyQt5.QtCore import Qt


class SNPMultiPredictApp(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("PyrusGS-5K Engine")
        self.resize(900, 650)

        # {trait_name: model_package}
        self.models = {}
        self.last_result_df = None

        # ======= 总体布局 =======
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout()
        central.setLayout(main_layout)

        # 左侧控制区
        control_frame = QFrame()
        control_frame.setFrameShape(QFrame.StyledPanel)
        control_layout = QVBoxLayout()
        control_frame.setLayout(control_layout)

        # 右侧日志区
        right_frame = QFrame()
        right_frame.setFrameShape(QFrame.StyledPanel)
        right_layout = QVBoxLayout()
        right_frame.setLayout(right_layout)

        main_layout.addWidget(control_frame, 4)
        main_layout.addWidget(right_frame, 6)

        # ---------- 1. 模型目录 ----------
        lbl_models_dir = QLabel("Model File")
        self.edit_models_dir = QLineEdit()
        # TODO：改成你自己的模型目录
        self.edit_models_dir.setText(r"D:/vir-nir/genes-type/pre_model")
        btn_browse_dir = QPushButton("Files")
        btn_browse_dir.clicked.connect(self.browse_models_dir)
        btn_load_all = QPushButton("Load all models")
        btn_load_all.clicked.connect(self.load_all_models)

        row_dir_1 = QHBoxLayout()
        row_dir_1.addWidget(self.edit_models_dir)
        row_dir_1.addWidget(btn_browse_dir)

        control_layout.addWidget(lbl_models_dir)
        control_layout.addLayout(row_dir_1)
        control_layout.addWidget(btn_load_all)

        # ---------- 2. 性状多选列表 ----------
        lbl_traits = QLabel("Select the traits to be predicted:")
        self.list_traits = QListWidget()
        self.list_traits.setSelectionMode(QListWidget.ExtendedSelection)
        control_layout.addWidget(lbl_traits)
        control_layout.addWidget(self.list_traits)

        # ---------- 3. 待预测 SNP CSV ----------
        lbl_data = QLabel("CSV file of SNPs to be predicted:")
        self.edit_data = QLineEdit()
        btn_browse_data = QPushButton("Select CSV")
        btn_browse_data.clicked.connect(self.browse_data)

        row_data = QHBoxLayout()
        row_data.addWidget(self.edit_data)
        row_data.addWidget(btn_browse_data)

        control_layout.addWidget(lbl_data)
        control_layout.addLayout(row_data)

        # ---------- 4. 按钮 ----------
        btn_row = QHBoxLayout()
        self.btn_predict = QPushButton("Start prediction (batch output)")
        self.btn_predict.clicked.connect(self.run_predict)
        self.btn_save = QPushButton("Save result as CSV")
        self.btn_save.clicked.connect(self.save_results)

        btn_row.addWidget(self.btn_predict)
        btn_row.addWidget(self.btn_save)
        control_layout.addLayout(btn_row)

        control_layout.addStretch(1)

        # ---------- 5. 右侧：预测结果文本框 ----------
        self.result_box = QTextEdit()
        self.result_box.setReadOnly(True)
        right_layout.addWidget(self.result_box)

        # ---------- 6. 右侧：日志 ----------
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setFixedHeight(200)
        right_layout.addWidget(self.log_box, 3)

        # ---------- 7. 样式 ----------
        self.apply_style()

    # ========== 美化样式 ==========
    def apply_style(self):
        qss = """
        QWidget {
            font-family: 'Microsoft YaHei';
            font-size: 11pt;
        }
        QMainWindow {
            background-color: #f5f7fb;
        }
        QFrame {
            background-color: #ffffff;
            border-radius: 10px;
        }
        QLineEdit, QTextEdit, QListWidget {
            border: 1px solid #d0d4e6;
            border-radius: 6px;
            padding: 4px;
            background-color: #fbfcff;
        }
        QPushButton {
            border-radius: 6px;
            padding: 6px 14px;
            border: 1px solid #4c6fff;
            background-color: #4c6fff;
            color: #ffffff;
        }
        QPushButton:hover {
            background-color: #3d59d6;
        }
        QPushButton:pressed {
            background-color: #3347b3;
        }
        QLabel {
            color: #333333;
        }
        """
        self.setStyleSheet(qss)

    # ========== 日志 ==========
    def log(self, msg: str):
        self.log_box.append(msg)
        self.log_box.moveCursor(self.log_box.textCursor().End)
        print(msg)

    # ========== 选择模型目录 ==========
    def browse_models_dir(self):
        d = QFileDialog.getExistingDirectory(self, "选择模型目录")
        if d:
            self.edit_models_dir.setText(d)

    # ========== 加载所有模型 ==========
    def load_all_models(self):
        models_dir = self.edit_models_dir.text().strip()
        if not models_dir or not os.path.isdir(models_dir):
            QMessageBox.warning(self, "提示", "请先选择正确的模型目录。")
            return

        pattern = os.path.join(models_dir, "new_final_model_*.pt")
        files = glob.glob(pattern)
        if not files:
            QMessageBox.warning(self, "提示", f"目录中未找到 final_model_*.pt：\n{models_dir}")
            return

        self.models.clear()
        self.list_traits.clear()
        self.log_box.clear()

        self.log(f"在目录中找到 {len(files)} 个模型文件，开始加载...")

        loaded = 0
        for path in files:
            fname = os.path.basename(path)
            try:
                with open(path, "rb") as f:
                    pkg = pickle.load(f)

                if "stacking_model" not in pkg or "feature_columns" not in pkg:
                    self.log(f"[跳过] {fname} 缺少 stacking_model 或 feature_columns")
                    continue

                trait = pkg.get("target_column", None)
                if not trait:
                    trait = fname.replace("final_model_", "").replace(".pt", "")

                self.models[trait] = pkg
                loaded += 1
                self.log(f"[已加载] {fname}  -> 性状：{trait}")
            except Exception as e:
                self.log(f"[失败] 加载 {fname} 时出错：{e}")

        if loaded == 0:
            QMessageBox.critical(self, "错误", "没有成功加载任何模型，请检查 pt 文件。")
            return

        # 多选列表（带勾选框）
        for trait in sorted(self.models.keys()):
            item = QListWidgetItem(trait)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked)
            self.list_traits.addItem(item)

        self.log(f"\n共成功加载 {loaded} 个性状模型。")
        QMessageBox.information(self, "完成", f"成功加载 {loaded} 个模型。")

    # ========== 选择 SNP CSV ==========
    def browse_data(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "选择待预测 SNP 数据 CSV", "",
            "CSV Files (*.csv);;All Files (*)"
        )
        if file_path:
            self.edit_data.setText(file_path)

    # ========== 批量预测并直接显示预测结果 ==========
    def run_predict(self):
        if not self.models:
            QMessageBox.warning(self, "提示", "请先加载模型目录。")
            return

        # 勾选的性状
        selected_traits = []
        for i in range(self.list_traits.count()):
            item = self.list_traits.item(i)
            if item.checkState() == Qt.Checked:
                selected_traits.append(item.text())

        if not selected_traits:
            QMessageBox.warning(self, "提示", "请至少勾选一个需要预测的性状。")
            return

        data_path = self.edit_data.text().strip()
        if not data_path:
            QMessageBox.warning(self, "提示", "请先选择待预测 SNP CSV 文件。")
            return

        # 读数据
        try:
            df = pd.read_csv(data_path)
        except Exception as e:
            self.log(f"[错误] 读取 CSV 失败：{e}")
            QMessageBox.critical(self, "错误", f"读取 CSV 失败：\n{e}")
            return

        df_result = df.copy()
        self.result_box.clear()
        self.log_box.clear()
        self.log(f"开始对 {len(df)} 条样本进行多性状预测...")
        self.log(f"选择的性状：{', '.join(selected_traits)}\n")

        result_str = ""

        # 逐个性状预测并实时显示结果
        for trait in selected_traits:
            pkg = self.models.get(trait)
            if pkg is None:
                self.log(f"[跳过] 未找到性状 {trait} 的模型")
                continue

            if "pca" not in pkg:
                self.log(f"[错误] {trait} 模型未保存 'pca' 对象，请在训练脚本中一并保存。")
                continue

            scaler = pkg.get("scaler", None)
            pca = pkg["pca"]
            selected_features = pkg.get("selected_features", None)
            model = pkg["stacking_model"]
            feature_cols = pkg["feature_columns"]

            # 对齐特征列
            missing = [c for c in feature_cols if c not in df.columns]
            if missing:
                self.log(f"[跳过] {trait}：数据缺少特征列 {missing}")
                continue

            X = df[feature_cols].values

            try:
                # 标准化
                if scaler is not None:
                    X_scaled = scaler.transform(X)
                else:
                    X_scaled = X

                # PCA
                X_pca = pca.transform(X_scaled)

                # 特征选择
                if selected_features is not None:
                    X_sel = X_pca[:, selected_features]
                else:
                    X_sel = X_pca

                # 预测（分类/回归都兼容）
                y_pred = model.predict(X_sel)

            except Exception as e:
                self.log(f"[错误] {trait} 预测失败：{e}")
                continue

            pred_col = f"Predicted_{trait}"
            df_result[pred_col] = y_pred
            result_str += f"{trait}: {y_pred[0]}\n"  # 这里只显示第一个样本的预测值
            self.log(f"[完成] {trait} 预测，结果列：{pred_col}")

        self.last_result_df = df_result
        self.log("\n全部性状预测已完成，将结果显示并导出为新的 CSV。")

        # 显示预测结果（按性状和预测值）
        self.result_box.append("预测结果：\n")
        self.result_box.append(result_str)

    # ========== 保存结果为 CSV ==========
    def save_results(self):
        if self.last_result_df is None:
            QMessageBox.information(self, "提示", "当前没有可保存的结果，请先进行预测。")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "保存预测结果为 CSV", "", "CSV Files (*.csv);;All Files (*)"
        )
        if not file_path:
            return

        try:
            self.last_result_df.to_csv(file_path, index=False)
            QMessageBox.information(self, "成功", f"预测结果已保存到：\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self, "错误", f"保存失败：\n{e}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SNPMultiPredictApp()
    window.show()
    sys.exit(app.exec_())
