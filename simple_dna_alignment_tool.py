import streamlit as st
import numpy as np

# 🎥 Video Background Styling
st.markdown(
    f"""
    <style>
    .video-background {{
        position: fixed;
        right: 0;
        bottom: 0;
        min-width: 100%;
        min-height: 100%;
        z-index: -1;
        object-fit: cover;
        opacity: 0.2;
    }}
    .block-container {{
        position: relative;
        z-index: 1;
        backdrop-filter: blur(4px); /* makes content stand out */
    }}
    </style>

    <video autoplay muted loop class="video-background">
        <source src="https://github.com/anushkagade/updatedproject/raw/main/background.mp4" type="video/mp4">
    </video>
    """,
    unsafe_allow_html=True
)

# 🧬 Scoring system
MATCH_AWARD = 1
MISMATCH_PENALTY = -1
GAP_PENALTY = -2

# 🔬 Needleman-Wunsch Algorithm
def needleman_wunsch(seq1, seq2):
    m, n = len(seq1), len(seq2)
    score = np.zeros((m + 1, n + 1), dtype=int)

    for i in range(m + 1):
        score[i][0] = GAP_PENALTY * i
    for j in range(n + 1):
        score[0][j] = GAP_PENALTY * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + (MATCH_AWARD if seq1[i - 1] == seq2[j - 1] else MISMATCH_PENALTY)
            delete = score[i - 1][j] + GAP_PENALTY
            insert = score[i][j - 1] + GAP_PENALTY
            score[i][j] = max(match, delete, insert)

    align1, align2 = '', ''
    i, j = m, n

    while i > 0 and j > 0:
        current = score[i][j]
        diagonal = score[i - 1][j - 1]
        up = score[i - 1][j]
        left = score[i][j - 1]

        if current == diagonal + (MATCH_AWARD if seq1[i - 1] == seq2[j - 1] else MISMATCH_PENALTY):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif current == up + GAP_PENALTY:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    while i > 0:
        align1 = seq1[i - 1] + align1
        align2 = '-' + align2
        i -= 1

    while j > 0:
        align1 = '-' + align1
        align2 = seq2[j - 1] + align2
        j -= 1

    return align1, align2, score[m][n]

# 🌐 Streamlit UI
st.title("🔬 Simple DNA Sequence Alignment Tool")

method = st.radio("Choose alignment method:", ["Needleman-Wunsch (Global)"])

seq1 = st.text_input("Enter DNA Sequence 1 (e.g. AGCTG):")
seq2 = st.text_input("Enter DNA Sequence 2 (e.g. AGCT):")

if st.button("Align Sequences"):
    if not seq1 or not seq2:
        st.warning("Please enter both DNA sequences!")
    else:
        align1, align2, score = needleman_wunsch(seq1.upper(), seq2.upper())
        st.subheader("Alignment Result")
        st.text(f"Sequence 1: {align1}")
        st.text(f"Sequence 2: {align2}")
        st.success(f"Alignment Score: {score}")
