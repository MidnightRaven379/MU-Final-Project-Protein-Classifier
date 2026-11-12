from pathlib import Path
from urllib.request import Request, urlopen


# ==========================================
# DOWNLOAD PROTEOME FROM UNIPROT
# ==========================================

def download_protein_data():

    url = (
        "https://rest.uniprot.org/uniprotkb/stream?"
        "compressed=false"
        "&format=tsv"
        "&query=proteome:UP000005640"
        "&fields=accession,protein_name,protein_families,sequence"
    )

    try:

        request = Request(
            url,
            headers={"User-Agent": "Python"}
        )

        with urlopen(request) as response:
            data = response.read().decode("utf-8")

        proteins = {}

        lines = data.splitlines()

        for line in lines[1:]:

            fields = line.split("\t")

            if len(fields) < 4:
                continue

            accession = fields[0]
            protein_name = fields[1]
            protein_class = fields[2]
            sequence = fields[3]

            if protein_class.strip() == "":
                protein_class = "Unclassified"

            proteins[accession] = {
                "name": protein_name,
                "class": protein_class,
                "sequence": sequence
            }

        return proteins

    except Exception as e:

        print("Download failed:")
        print(e)

        return {}


# ==========================================
# RNA -> PROTEIN TRANSLATION
# ==========================================

def translate_rna(rna):

    rna = rna.upper().replace("T", "U")

    valid_bases = {"A", "U", "G", "C"}

    for base in rna:
        if base not in valid_bases:
            return None

    codon_table = {

        "UUU": "F", "UUC": "F",
        "UUA": "L", "UUG": "L",

        "CUU": "L", "CUC": "L",
        "CUA": "L", "CUG": "L",

        "AUU": "I", "AUC": "I",
        "AUA": "I",

        "AUG": "M",

        "GUU": "V", "GUC": "V",
        "GUA": "V", "GUG": "V",

        "UCU": "S", "UCC": "S",
        "UCA": "S", "UCG": "S",

        "CCU": "P", "CCC": "P",
        "CCA": "P", "CCG": "P",

        "ACU": "T", "ACC": "T",
        "ACA": "T", "ACG": "T",

        "GCU": "A", "GCC": "A",
        "GCA": "A", "GCG": "A",

        "UAU": "Y", "UAC": "Y",

        "CAU": "H", "CAC": "H",

        "CAA": "Q", "CAG": "Q",

        "AAU": "N", "AAC": "N",

        "AAA": "K", "AAG": "K",

        "GAU": "D", "GAC": "D",

        "GAA": "E", "GAG": "E",

        "UGU": "C", "UGC": "C",

        "UGG": "W",

        "CGU": "R", "CGC": "R",
        "CGA": "R", "CGG": "R",

        "AGU": "S", "AGC": "S",

        "AGA": "R", "AGG": "R",

        "GGU": "G", "GGC": "G",
        "GGA": "G", "GGG": "G",

        "UAA": "*",
        "UAG": "*",
        "UGA": "*"
    }

    protein = ""

    for i in range(0, len(rna) - 2, 3):

        codon = rna[i:i+3]

        amino_acid = codon_table.get(codon)

        if amino_acid is None:
            continue

        if amino_acid == "*":
            break

        protein += amino_acid

    return protein


# ==========================================
# SEARCH BY NAME
# ==========================================

def search_by_name(proteins):

    keyword = input(
        "\nEnter Protein Name: "
    ).lower().strip()

    results = []

    for accession, info in proteins.items():

        if keyword in info["name"].lower():

            results.append(
                (accession, info)
            )

    if len(results) == 0:

        print("\nNo proteins found.")
        return

    print(
        "\nFound",
        len(results),
        "matching proteins\n"
    )

    for accession, info in results[:10]:

        print("UniProt ID:",
              accession)

        print("Name:",
              info["name"])

        print("Class:",
              info["class"])

        print("Length:",
              len(info["sequence"]))

        print("-" * 60)


# ==========================================
# SEARCH BY RNA
# ==========================================

def search_by_rna(proteins):

    rna = input(
        "\nEnter RNA/DNA sequence: "
    ).strip()

    protein_fragment = translate_rna(rna)

    if protein_fragment is None:

        print("Invalid sequence entered.")
        return

    if len(protein_fragment) == 0:

        print("Translation produced no protein.")
        return

    matches = []

    for accession, info in proteins.items():

        if protein_fragment in info["sequence"]:

            matches.append(
                (accession, info)
            )

    if len(matches) == 0:

        print(
            "\nNo matching proteins found."
        )

        return

    print(
        "\nFound",
        len(matches),
        "matching proteins\n"
    )

    for accession, info in matches[:10]:

        print("UniProt ID:",
              accession)

        print("Name:",
              info["name"])

        print("Class:",
              info["class"])

        print("Length:",
              len(info["sequence"]))

        print("-" * 60)


# ==========================================
# LIST BY PROTEIN CLASS
# ==========================================

def search_by_class(proteins):

    keyword = input(
        "\nEnter protein class: "
    ).lower().strip()

    found = 0

    print("\nMatching Proteins:\n")

    for accession, info in proteins.items():

        if keyword in info["class"].lower():

            print(
                accession,
                "-",
                info["name"]
            )

            found += 1

    print(
        "\nTotal proteins found:",
        found
    )


# ==========================================
# DATABASE STATISTICS
# ==========================================

def database_statistics(proteins):

    total = len(proteins)

    if total == 0:

        print("Database empty.")
        return

    longest_name = ""
    longest_length = 0

    for info in proteins.values():

        length = len(info["sequence"])

        if length > longest_length:

            longest_length = length
            longest_name = info["name"]

    print("\nDATABASE STATISTICS")

    print("---------------------")

    print("Total Proteins:",
          total)

    print("Longest Protein:",
          longest_name)

    print("Longest Length:",
          longest_length)


# ==========================================
# SAVE SUMMARY
# ==========================================

def save_summary(proteins):

    Path("output").mkdir(
        exist_ok=True
    )

    file_path = (
        "output/protein_summary.txt"
    )

    with open(
        file_path,
        "w",
        encoding="utf-8"
    ) as file:

        file.write(
            "Human Proteome Database\n"
        )

        file.write(
            f"Total Proteins: {len(proteins)}\n"
        )

    print(
        "\nSummary saved:"
    )

    print(file_path)


# ==========================================
# MAIN
# ==========================================

def main():

    print(
        "Downloading Human Proteome..."
    )

    proteins = download_protein_data()

    if len(proteins) == 0:

        print(
            "Failed to load data."
        )
        return

    print(
        "\nLoaded",
        len(proteins),
        "proteins."
    )

    while True:

        print(
            "\n===== HUMAN PROTEOME DATABASE ====="
        )

        print(
            "1. Search Protein by Name"
        )

        print(
            "2. Search Protein by RNA Sequence"
        )

        print(
            "3. Search Proteins by Class"
        )

        print(
            "4. Database Statistics"
        )

        print(
            "5. Save Summary"
        )

        print(
            "6. Exit"
        )

        choice = input(
            "\nEnter Choice: "
        )

        if choice == "1":

            search_by_name(
                proteins
            )

        elif choice == "2":

            search_by_rna(
                proteins
            )

        elif choice == "3":

            search_by_class(
                proteins
            )

        elif choice == "4":

            database_statistics(
                proteins
            )

        elif choice == "5":

            save_summary(
                proteins
            )

        elif choice == "6":

            print("Goodbye!")
            break

        else:

            print(
                "Invalid choice."
            )


if __name__ == "__main__":
    main()
