import streamlit as st
from protein_classifier import (
    download_protein_data,
    translate_rna
)

# ----------------------------------
# LOAD DATA
# ----------------------------------
@st.cache_data
def load_proteins():
    return download_protein_data()

proteins = load_proteins()

# ----------------------------------
# PAGE SETTINGS
# ----------------------------------
st.set_page_config(
    page_title="Protein Classifier",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Human Protein Classifier")

st.write(
    f"Database loaded successfully with {len(proteins)} proteins."
)

menu = st.sidebar.selectbox(
    "Choose an Option",
    [
        "Search by Name",
        "Search by RNA",
        "Search by Class",
        "Database Statistics",
        "Save Summary"
    ]
)

# ----------------------------------
# SEARCH BY NAME
# ----------------------------------
if menu == "Search by Name":

    st.header("Search Protein by Name")

    keyword = st.text_input(
        "Enter Protein Name"
    )

    if st.button("Search Name"):

        results = []

        for accession, info in proteins.items():

            if keyword.lower() in info["name"].lower():

                results.append(
                    (accession, info)
                )

        st.success(
            f"{len(results)} proteins found"
        )

        for accession, info in results[:10]:

            st.write(f"**UniProt ID:** {accession}")
            st.write(f"**Name:** {info['name']}")
            st.write(f"**Class:** {info['class']}")
            st.write(
                f"**Length:** {len(info['sequence'])}"
            )

            st.divider()

# ----------------------------------
# SEARCH BY RNA
# ----------------------------------
elif menu == "Search by RNA":

    st.header("Search Protein by RNA")

    rna = st.text_area(
        "Enter RNA/DNA Sequence"
    )

    if st.button("Search RNA"):

        fragment = translate_rna(rna)

        if fragment is None:

            st.error(
                "Invalid RNA/DNA sequence."
            )

        else:

            matches = []

            for accession, info in proteins.items():

                if fragment in info["sequence"]:

                    matches.append(
                        (accession, info)
                    )

            st.success(
                f"{len(matches)} matches found"
            )

            for accession, info in matches[:10]:

                st.write(
                    f"**UniProt ID:** {accession}"
                )
                st.write(
                    f"**Name:** {info['name']}"
                )
                st.write(
                    f"**Class:** {info['class']}"
                )

                st.divider()

# ----------------------------------
# SEARCH BY CLASS
# ----------------------------------
elif menu == "Search by Class":

    st.header("Search Protein by Class")

    keyword = st.text_input(
        "Enter Protein Class"
    )

    if st.button("Search Class"):

        found = []

        for accession, info in proteins.items():

            if keyword.lower() in info["class"].lower():

                found.append(
                    (accession, info)
                )

        st.success(
            f"{len(found)} proteins found"
        )

        for accession, info in found[:50]:

            st.write(
                f"{accession} | {info['name']}"
            )

# ----------------------------------
# DATABASE STATS
# ----------------------------------
elif menu == "Database Statistics":

    st.header("Database Statistics")

    total = len(proteins)

    longest_name = ""
    longest_length = 0

    for info in proteins.values():

        length = len(info["sequence"])

        if length > longest_length:

            longest_length = length
            longest_name = info["name"]

    st.metric(
        "Total Proteins",
        total
    )

    st.metric(
        "Longest Protein Length",
        longest_length
    )

    st.write(
        "Longest Protein:",
        longest_name
    )

# ----------------------------------
# SAVE SUMMARY
# ----------------------------------
elif menu == "Save Summary":

    st.header("Save Summary")

    if st.button("Generate Summary"):

        with open(
            "protein_summary.txt",
            "w",
            encoding="utf-8"
        ) as file:

            file.write(
                "Human Proteome Database\n"
            )

            file.write(
                f"Total Proteins: {len(proteins)}\n"
            )

        st.success(
            "Summary saved as protein_summary.txt"
        )
