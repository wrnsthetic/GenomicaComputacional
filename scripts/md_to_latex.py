"""Converte relatorio/relatorio.md em um projeto LaTeX pronto para o Overleaf (relatorio/latex/).

Uso: python scripts/md_to_latex.py   (precisa do pandoc no PATH)
Saída: relatorio/latex/main.tex + relatorio/latex/figuras/*.png  → compactar a pasta e subir no Overleaf (pdfLaTeX).
"""
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MD = ROOT / "relatorio" / "relatorio.md"
OUT = ROOT / "relatorio" / "latex"
FIG_SRC = ROOT / "results" / "figures"

# símbolos que o pdfLaTeX não aceita direto no texto
SIMBOLOS = {"≥": r"$\geq$", "≤": r"$\leq$", "≈": r"$\approx$", "≫": r"$\gg$", "↔": r"$\leftrightarrow$",
            "→": r"$\rightarrow$", "×": r"$\times$", "−": r"$-$", "√": r"$\surd$", "∈": r"$\in$", "⚠️": "", "·": r"$\cdot$"}

PREAMBULO = r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage[margin=2.5cm]{geometry}
\usepackage{graphicx}
\usepackage{booktabs,longtable,array,calc}
\usepackage{float}
\usepackage{amsmath,amssymb}
\usepackage[hidelinks]{hyperref}
\usepackage{microtype}
\setlength{\parskip}{0.5em}
\setlength{\parindent}{0pt}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\graphicspath{{figuras/}}

\title{%TITULO%}
\author{Laís Goulart\\[2pt]\small Disciplina: Genômica Computacional}
\date{\today}

\begin{document}
\maketitle

\noindent\textbf{Artigo de referência:} Liu, M. C. et al. (2020). \emph{Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA}. Annals of Oncology, 31(6), 745--759.\\
\textbf{Notebook reproduzível:} \texttt{notebooks/trabalho\_final.ipynb}

"""


def figura(m):
    arquivo, legenda = m.group(1).strip(), m.group(2).strip()
    return ("\n```{=latex}\n\\begin{figure}[H]\n\\centering\n"
            f"\\includegraphics[width=0.85\\linewidth]{{{arquivo}}}\n\\caption{{{legenda}}}\n\\end{{figure}}\n```\n")


def preparar_md(t):
    titulo = re.search(r"^# (.+)$", t, re.M).group(1)
    t = t[t.index("## Resumo"):]                                   # descarta título e cabeçalho
    t = re.sub(r"^(#{2,}) \d+(?:\.\d+)*\.? ", r"\1 ", t, flags=re.M)  # numeração fica com o LaTeX
    t = re.sub(r"^## (Registro das fontes.*|Referências)$", r"## \1 {-}", t, flags=re.M)
    t = re.sub(r"^\*\*(\d+\.\d+) ", "**", t, flags=re.M)            # "**6.1 Pergunta**" → "**Pergunta**"
    t = re.sub(r"`\[Figura: ([^|]+)\|([^\]]+)\]`", figura, t)
    t = t.replace("\n---\n", "\n")
    # pandoc só reconhece lista depois de linha em branco
    t = re.sub(r"(?m)^((?!\s*(?:[-*] |\d+\. |\|)).+)\n(?=(?:- |\d+\. ))", r"\1\n\n", t)
    # resumo → environment abstract
    a, b = t.index("## Resumo"), t.index("## Introdução")
    resumo = t[a + len("## Resumo"):b].strip()
    t = "```{=latex}\n\\begin{abstract}\n```\n" + resumo + "\n```{=latex}\n\\end{abstract}\n```\n\n" + t[b:]
    return titulo, t


def main():
    titulo, corpo_md = preparar_md(MD.read_text(encoding="utf-8"))
    corpo = subprocess.run(["pandoc", "-f", "markdown-auto_identifiers", "-t", "latex", "--columns=80",
                            "--shift-heading-level-by=-1"],
                           input=corpo_md, capture_output=True, text=True, encoding="utf-8", check=True).stdout
    for s, r in SIMBOLOS.items():
        corpo = corpo.replace(s, r)
    estranhos = sorted({c for c in corpo if ord(c) > 255} - {"–", "—"})  # travessões o inputenc aceita
    if estranhos:
        print("ATENCAO, caracteres fora do Latin-1:", [hex(ord(c)) for c in estranhos])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figuras").mkdir(exist_ok=True)
    for f in re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", corpo):
        shutil.copy(FIG_SRC / f, OUT / "figuras" / f)
    (OUT / "main.tex").write_text(PREAMBULO.replace("%TITULO%", titulo) + corpo + "\n\\end{document}\n", encoding="utf-8")
    print("escrito:", OUT / "main.tex")


if __name__ == "__main__":
    main()
