# ==================================================
# File: rules/0_prepare_ogm_data.smk
# Description:
#   Prepares Optical Genome Mapping (OGM) input data for downstream analysis.
#   Depending on the configuration (ZIP = True/False), this rule either:
#     1️⃣ Waits for decompressed outputs produced by `decompressOGM` (per patient)
#     2️⃣ Or copies raw OGM data directly from the input.bionanodata directory
#
# Author: Marta Portasany
# Updated: 2025-11-05
# Pipeline: PULPO 🐙
# Requirements:
#   - The following variables must be defined in the main Snakefile:
#       config, sample_table, resultsdir, logsdir
# ==================================================

ZIP = bool(config.get('analysis', {}).get('zip', False))

if ZIP:
    # Case ZIP=True: decompression is handled by `decompressOGM` (per patient),
    # which generates .../OGMdata/.done.
    # Here, we simply aggregate those completion markers and create our global done.txt.
    rule prepare_ogm_data:
        input:
            descomp = expand(
                f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata/.done",
                anonymised=sample_table["anonymised"]
            )
        output:
            done = f"{logsdir}/DATA/prepare_ogm_data/done.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.done})"
            echo ok > "{output.done}"
            """

else:
    # Case ZIP=False: copy only relevant OGM files (.smap and any file containing 'CNV')
    # from input.bionanodata/<sample> to results/DATA/Patients/<anonymised>/OGMdata,
    # creating all required directories.
    rule prepare_ogm_data:
        input:
            # no input dependencies
        output:
            done = f"{logsdir}/DATA/prepare_ogm_data/done.txt",
            ogm_dir = directory(
                expand(
                    f"{resultsdir}/DATA/Patients/{{anonymised}}/OGMdata",
                    anonymised=sample_table["anonymised"]
                )
            )
        run:
            import os, shutil, pathlib, re

            base_src = config["input"]["bionanodata"]

            smap_pattern = re.compile(r"\.smap$")
            cnv_pattern  = re.compile(r"CNV", re.IGNORECASE)

            for _, row in sample_table.iterrows():
                anon = row["anonymised"]
                sample_id = row["sample"]

                source = os.path.join(base_src, sample_id)
                target = f"{resultsdir}/DATA/Patients/{anon}/OGMdata"

                pathlib.Path(target).mkdir(parents=True, exist_ok=True)

                if not os.path.exists(source):
                    raise FileNotFoundError(f"OGM data not found: {source}")

                copied_files = 0

                for fname in os.listdir(source):
                    src_f = os.path.join(source, fname)
                    dst_f = os.path.join(target, fname)

                    # Keep only files/folders that match .smap or contain 'CNV'
                    if smap_pattern.search(fname) or cnv_pattern.search(fname):
                        if os.path.isdir(src_f):
                            if os.path.exists(dst_f):
                                shutil.rmtree(dst_f)
                            shutil.copytree(src_f, dst_f)
                        else:
                            shutil.copy2(src_f, dst_f)
                        copied_files += 1

                print(f"✅ Copied {copied_files} files for {anon} ({sample_id})")

            # Write global 'done' marker
            pathlib.Path(os.path.dirname(output.done)).mkdir(parents=True, exist_ok=True)
            with open(output.done, "w") as fh:
                fh.write("ok\n")

