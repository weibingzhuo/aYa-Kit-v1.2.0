import tkinter as tk
from tkinter import ttk, scrolledtext
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def calculate_similarity(alignment):
    seqA = alignment.seqA.replace("", "")
    seqB = alignment.seqB.replace("-", "")
    matches = sum(a == b for a, b in zip(seqA, seqB))
    total = max(len(seqA), len(seqB))
    return matches / total * 100 if total != 0 else 0


def kabat_rule(seq):
    # 示例Kabat规则（假设的区域切分）
    return {
        "FR1": (0, 23),
        "CDR1": (23, 31),
        "FR2": (31, 48),
        "CDR2": (48, 56),
        "FR3": (56, 71),
        "CDR3": (71, 88),
        "FR4": (88, len(seq))
    }


def imgt_rule(seq):
    # 示例IMGT规则（假设的区域切分）
    return {
        "FR1": (0, 26),
        "CDR1": (26, 34),
        "FR2": (34, 54),
        "CDR2": (54, 61),
        "FR3": (61, 85),
        "CDR3": (85, 100),
        "FR4": (100, len(seq))
    }


def chothia_rule(seq):
    # 示例Chothia规则（假设的区域切分）
    return {
        "FR1": (0, 24),
        "CDR1": (24, 32),
        "FR2": (32, 53),
        "CDR2": (53, 61),
        "FR3": (61, 82),
        "CDR3": (82, 98),
        "FR4": (98, len(seq))
    }


def honneger_rule(seq):
    # 示例Honneger规则（假设的区域切分）
    return {
        "FR1": (0, 25),
        "CDR1": (25, 33),
        "FR2": (33, 55),
        "CDR2": (55, 65),
        "FR3": (65, 75),
        "CDR3": (75, 87),
        "FR4": (87, len(seq))
    }


class SequenceApp:
    def __init__(self, rot):
        self.root = rot
        self.root.title("aYa-Kit v1.2.0           Author：Wey_Benz")
        self.root.geometry("900x600")

        # 输入区域
        self.input_frame = ttk.Frame(self.root)
        self.input_frame.pack(fill=tk.X, padx=10, pady=5)

        # 添加初始输入框
        self.entries = []
        self.add_sequence_pair()  # 第一组输入框
        self.add_sequence_pair()  # 第二组输入框

        # 控制按钮与选项
        self.btn_frame = ttk.Frame(self.root)
        self.btn_frame.pack(pady=5)

        ttk.Button(self.btn_frame, text="Comparing", command=self.run_alignment).pack(side=tk.LEFT, padx=5)
        ttk.Button(self.btn_frame, text="CDR/FR", command=self.annotate_cdr_fr).pack(side=tk.LEFT, padx=5)

        # 选项菜单
        self.annotation_rule = tk.StringVar(value="Kabat")
        self.annotation_menu = ttk.Combobox(self.btn_frame, textvariable=self.annotation_rule, state="readonly")
        self.annotation_menu['values'] = ["Kabat", "IMGT", "Chothia", "Honneger"]
        self.annotation_menu.pack(side=tk.LEFT, padx=5)

        # 结果展示区域
        self.result_frame = ttk.Frame(self.root)
        self.result_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.alignment_result = scrolledtext.ScrolledText(self.result_frame, wrap=tk.NONE, height=15)
        self.alignment_result.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.annotation_result = scrolledtext.ScrolledText(self.result_frame, wrap=tk.NONE, height=15)
        self.annotation_result.pack(fill=tk.BOTH, expand=True)

        # 配置文本颜色
        self.alignment_result.tag_configure("gap", foreground="red")
        self.alignment_result.tag_configure("match", foreground="blue")
        self.annotation_result.tag_configure("cdr", background="yellow")
        self.annotation_result.tag_configure("fr", background="light blue")

    def add_sequence_pair(self):
        frame = ttk.Frame(self.input_frame)
        frame.pack(fill=tk.X, pady=2)

        ttk.Label(frame, text="Name:").pack(side=tk.LEFT)
        name_entry = ttk.Entry(frame, width=15)
        name_entry.pack(side=tk.LEFT, padx=5)

        ttk.Label(frame, text="Sequence:").pack(side=tk.LEFT)
        seq_entry = ttk.Entry(frame, width=80)
        seq_entry.pack(side=tk.LEFT, padx=5)

        self.entries.append((name_entry, seq_entry))

    def run_alignment(self):
        sequences = []
        for name_entry, seq_entry in self.entries:
            name = name_entry.get()
            seq = seq_entry.get().upper()
            if name and seq:
                sequences.append((name, seq))

        if len(sequences) < 2:
            self.alignment_result.insert(tk.END, "Please enter at least two sequences.\n")
            return

        self.alignment_result.delete(1.0, tk.END)

        alignments = pairwise2.align.localxx(sequences[0][1], sequences[1][1], score_only=False)
        if alignments:
            alignment = alignments[0]
            formatted = format_alignment(*alignment)
            self.colorize_alignment(formatted)
        else:
            self.alignment_result.insert(tk.END, "No valid comparison results were found.\n")

    def colorize_alignment(self, text):
        self.alignment_result.delete(1.0, tk.END)
        for line in text.split('\n'):
            if not line.strip():
                self.alignment_result.insert(tk.END, line + "\n")
                continue
            for char in line:
                if char == "-":
                    self.alignment_result.insert(tk.END, char, "gap")
                elif char.isupper() and line.strip()[0] != "|" and char != "|":
                    self.alignment_result.insert(tk.END, char, "match")
                else:
                    self.alignment_result.insert(tk.END, char)
            self.alignment_result.insert(tk.END, "\n")

    def annotate_cdr_fr(self):
        selected_rule = self.annotation_rule.get()
        annotation = {
            "Kabat": kabat_rule,
            "IMGT": imgt_rule,
            "Chothia": chothia_rule,
            "Honneger": honneger_rule
        }

        sequences = []
        for name_entry, seq_entry in self.entries:
            name = name_entry.get()
            seq = seq_entry.get().upper()
            if name and seq:
                sequences.append((name, seq))

        if not sequences:
            self.annotation_result.insert(tk.END, "Please enter the antibody sequence.\n")
            return

        self.annotation_result.delete(1.0, tk.END)

        for name, seq in sequences:
            self.annotation_result.insert(tk.END, f"{name} sequence:\n", "bold")
            regions = annotation[selected_rule](seq)
            self.display_regions(seq, regions, name)

    def display_regions(self, seq, regions, name):
        self.annotation_result.insert(tk.END, f"{name}\n", "bold")
        tags = []

        # 添加标签信息
        for region, (start, end) in regions.items():
            if region.startswith("CDR"):
                tags.append((start, end, "cdr"))
            else:
                tags.append((start, end, "fr"))

        # 标记序列
        current_char = 0
        for tag in sorted(tags, key=lambda x: x[0]):
            start, end, label = tag
            if start > current_char:
                self.annotation_result.insert(tk.END, seq[current_char:start], "normal")
            self.annotation_result.insert(tk.END, seq[start:end], label)
            current_char = end

        self.annotation_result.insert(tk.END, "\n\n", "normal")


if __name__ == "__main__":
    root = tk.Tk()
    app = SequenceApp(root)
    root.mainloop()
