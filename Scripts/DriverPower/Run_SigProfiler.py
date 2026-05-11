from SigProfilerAssignment import Analyzer as Analyze

def main_function():
    data = "/project/home/p201120/ryan/SPA_nunes/vcfs"

    Analyze.cosmic_fit(
        samples=data,
        output="/project/home/p201120/ryan/SPA_nunes/output_SBS96",
        input_type="vcf",
        context_type="96",
        genome_build="GRCh38",
        cosmic_version=3.5,
        cpu=128,
        verbose=True,
        make_plots=False
    )

if __name__ == "__main__":
    main_function()

