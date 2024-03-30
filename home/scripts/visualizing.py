import matplotlib.pyplot as plt

# Read the BLAST output file and extract E-values
evalues = []
with open("blast.out", "r") as f:
    for line in f:
        if line.startswith("  E-value"):
            evalue = float(line.split()[1])
            evalues.append(evalue)

# Create a histogram of the E-values
plt.hist(evalues)
plt.xlabel("E-value")
plt.ylabel("Count")
plt.title("Distribution of E-values in BLAST results")
plt.show()

# You can create other visualizations based on your needs (e.g., scatter plots of alignment scores)
# ...
