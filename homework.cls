\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{homework}[2023/02/27 Custom class for my homework]

\LoadClass{article}
\author{}\date{}

\RequirePackage[a4paper]{geometry}
\RequirePackage{graphicx}

\RequirePackage[T1]{fontenc}
\RequirePackage{mlmodern}

\RequirePackage{amsmath, amssymb, amsthm}
\RequirePackage{etoolbox, enumitem}
\theoremstyle{definition}
\newtheorem{problem}{Problem}
\renewcommand{\qedsymbol}{\(\lozenge\)}
\AtEndEnvironment{problem}{\qed}

\RequirePackage{tikz}

\RequirePackage{lstfiracode}
\RequirePackage[framed]{matlab-prettifier}
\lstnewenvironment{matlab}{
    \lstset{
        style=Matlab-editor,
        escapechar=`
    }
}{}
\newcommand{\matlabfile}[1]{
    \lstinputlisting[
        style=Matlab-editor,
        basicstyle=\small\ttfamily,
        escapechar=`
    ]{#1}
}

\newcommand*{\Matlab}{\textsc{Matlab}}
% \d originally means "under dot"
\renewcommand*{\d}[1][\,]{\ensuremath{\mathop{\textnormal{d}}\!#1}}
