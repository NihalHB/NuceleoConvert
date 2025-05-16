# NucleoConvert 

**NucleoConvert** is an all-in-one R Shiny application designed to seamlessly convert DNA sequences into RNA while unlocking a suite of quality checks, statistical insights, and rich interactive visualizations. Whether you're a researcher, educator, or student, this tool brings sequence exploration to your fingertips—with style.

---

## Features

- **DNA-to-RNA Conversion**:
  - Convert DNA sequences to RNA with optional cleaning (remove gaps, N's, and invalid characters).
  - Preserve or clean sequences based on user preferences.

- **Quality Control**:
  - Sequence length validation.
  - GC content analysis with warnings for high/low GC.
  - Detection of gaps, homopolymers, and invalid characters.

- **Advanced RNA Analysis**:
  - **Codon Usage**: Analyze codon frequencies and generate interactive visualizations.
  - **ORF Finder**: Identify open reading frames (ORFs) in RNA sequences.
  - **Poly-A Tail Detection**: Detect and analyze poly-A tails in mRNA sequences.

- **Interactive Visualizations**:
  - GC content pie chart.
  - Nucleotide frequency bar chart.
  - Dinucleotide frequency heatmap.
  - Positional nucleotide distribution.
  - Sequence complexity analysis.
  - Homopolymer distribution.

- **Data Export**:
  - Download RNA sequences in FASTA format.
  - Export sequence statistics and QC reports in CSV format.

---

## Installation

### Prerequisites
- **R** (version 4.2.1 or higher)
- **RStudio** (recommended for ease of use)
- Required R packages: `shiny`, `shinythemes`, `stringr`, `plotly`, `ggplot2`, `DT`, `shinyFeedback`

### Steps
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/NucleoConvert-Analytics-Suite.git
   cd NucleoConvert-Analytics-Suite
   ```

2. Install the required R packages:
   ```R
   install.packages(c("shiny", "shinythemes", "stringr", "plotly", "ggplot2", "DT", "shinyFeedback"))
   ```

3. Run the Shiny app:
   ```R
   shiny::runApp()
   ```

---

## Usage

1. **Input DNA Sequence**:
   - Paste your DNA sequence into the input box.
   - Choose whether to auto-clean the sequence (remove gaps, N's, and invalid characters).

2. **Convert to RNA**:
   - Click the "Process Sequence" button to convert DNA to RNA.

3. **Explore Results**:
   - View the RNA sequence and detailed statistics in the **RNA Output** tab.
   - Use the **Analytics** tab to explore interactive visualizations.
   - Perform advanced RNA analysis (codon usage, ORF finding, poly-A tail detection) in the **RNA Analysis** tab.

4. **Download Results**:
   - Download the RNA sequence in FASTA format.
   - Export sequence statistics and QC reports in CSV format.

---

## Contributing

We welcome contributions! If you'd like to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. Commit your changes:
   ```bash
   git commit -m "Add your commit message here"
   ```
4. Push to your branch:
   ```bash
   git push origin feature/your-feature-name
   ```
5. Open a pull request and describe your changes.

---

## License

This project is licensed under the **MIT License**. 
See the [LICENSE](LICENSE) file for details.

---

## Contact

For questions, feedback, or collaboration opportunities, please contact:

- **Nihal Habib**  
  Email: nhabib@um6ss.ma  
  LinkedIn: (https://www.linkedin.com/in/nihal-habib-5546a4166/)


