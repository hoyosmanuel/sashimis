I. FAU
==================

.. code-block:: bash


  mkdir -p ~/Desktop/shashimis/gtf
  mkdir -p ~/Desktop/shashimis/svg
  
  python - <<'PY'
  from pathlib import Path
  import subprocess
  import sys
  
  ROOT = Path("/Volumes/Expansion/IGV")
  OUT = Path.home() / "Desktop" / "shashimis"
  GTF_DIR = OUT / "gtf"
  SVG_DIR = OUT / "svg"
  GGS = Path.home() / "herramental" / "ggsashimi" / "ggsashimi.py"
  
  GTF_DIR.mkdir(parents=True, exist_ok=True)
  SVG_DIR.mkdir(parents=True, exist_ok=True)
  
  def first_match(species, pattern):
      folder = ROOT / species
      matches = sorted(
          p for p in folder.glob(pattern)
          if p.is_file() and not p.name.startswith("._")
      )
      if len(matches) != 1:
          raise RuntimeError(
              f"Esperaba un archivo para {folder / pattern}; "
              f"encontré {len(matches)}: {[p.name for p in matches]}"
          )
      return matches[0]
  
  def grep_to_file(transcripts, infile, outfile):
      transcripts = set(transcripts)
      found = set()
  
      with infile.open() as source, outfile.open("w") as dest:
          for line in source:
              for transcript in transcripts:
                  if f'transcript_id "{transcript}"' in line:
                      dest.write(line)
                      found.add(transcript)
                      break
  
      missing = transcripts - found
      if missing:
          raise ValueError(
              f"Transcritos no encontrados en {infile}: {missing}"
          )
  
  def bed12_to_gtf(bed_path, out_path):
      records = []
  
      with bed_path.open() as f:
          for line in f:
              fields = line.rstrip("\n").split("\t")
              if len(fields) >= 12 and "#ACTB#" in fields[3]:
                  records.append(fields)
  
      if not records:
          raise ValueError(f"No encontré ACTB en {bed_path}")
  
      # Eliminar duplicados con estructura BED12 idéntica
      structures = {
          tuple(r[:3] + r[5:12]): r
          for r in records
      }
  
      if len(structures) != 1:
          raise ValueError(
              f"ACTB tiene {len(structures)} estructuras distintas "
              f"en {bed_path}. Requiere revisión."
          )
  
      fields = next(iter(structures.values()))
  
      chrom = fields[0]
      chrom_start = int(fields[1])
      chrom_end = int(fields[2])
      strand = fields[5]
      block_count = int(fields[9])
      block_sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
      block_starts = [int(x) for x in fields[11].rstrip(",").split(",")]
  
      if not (block_count == len(block_sizes) == len(block_starts)):
          raise ValueError(f"BED12 inconsistente: {bed_path}")
  
      with out_path.open("w") as out:
          attrs = 'gene_id "ACTB"; transcript_id "ACTB";'
          out.write(
              f"{chrom}\tTOGA\ttranscript\t"
              f"{chrom_start + 1}\t{chrom_end}\t.\t"
              f"{strand}\t.\t{attrs}\n"
          )
  
          for i, (size, offset) in enumerate(
              zip(block_sizes, block_starts), start=1
          ):
              start = chrom_start + offset + 1
              end = chrom_start + offset + size
  
              out.write(
                  f"{chrom}\tTOGA\texon\t"
                  f"{start}\t{end}\t.\t{strand}\t.\t"
                  f'gene_id "ACTB"; transcript_id "ACTB"; '
                  f'exon_number "{i}";\n'
              )
  
  def run_ggsashimi(bam, coords, gtf, out_prefix):
      cmd = [
          sys.executable, str(GGS),
          "-b", str(bam),
          "-c", coords,
          "-g", str(gtf),
          "-o", str(out_prefix),
          "-F", "svg",
          "--shrink",
          "-M", "1",
          "--width", "8",
          "--height", "2.5",
          "--ann-height", "1.4",
          "--base-size", "12"
      ]
  
      print(f"\nGenerando: {out_prefix.name}", flush=True)
      subprocess.run(cmd, check=True)
  
  # ---------------------------------
  # Archivos originales
  # ---------------------------------
  ESPE = "Eonycteris_spelaea"
  RAEG = "Rousettus_aegyptiacus"
  
  runs = {
      "SRR6958995": ESPE,
      "SRR6958999": ESPE,
      "SRR6959000": ESPE,
      "SRR18652228": RAEG
  }
  
  bams = {
      run: first_match(species, f"{run}*.bam")
      for run, species in runs.items()
  }
  
  gtfs = {
      run: first_match(species, f"{run}*.gtf")
      for run, species in runs.items()
  }
  
  bed_espe = first_match(ESPE, "*TOGA.bed")
  bed_raeg = first_match(RAEG, "*TOGA.bed")
  
  # ---------------------------------
  # GTF específicos de FAU
  # ---------------------------------
  fau_transcripts = {
      "SRR6958995": ["STRG.3514.1"],
      "SRR6958999": ["STRG.4849.2"],
      "SRR6959000": ["STRG.4297.1"],
      "SRR18652228": ["STRG.3989.1", "STRG.3989.3"]
  }
  
  fau_gtfs = {}
  
  for run, transcripts in fau_transcripts.items():
      species_prefix = "eSpe" if runs[run] == ESPE else "rAeg"
      destination = GTF_DIR / f"{species_prefix}_{run}_FAU.gtf"
  
      grep_to_file(transcripts, gtfs[run], destination)
      fau_gtfs[run] = destination
  
  # ---------------------------------
  # GTF de ACTB desde TOGA BED12
  # ---------------------------------
  actb_espe_gtf = GTF_DIR / "Eonycteris_spelaea_ACTB.gtf"
  actb_raeg_gtf = GTF_DIR / "Rousettus_aegyptiacus_ACTB.gtf"
  
  bed12_to_gtf(bed_espe, actb_espe_gtf)
  bed12_to_gtf(bed_raeg, actb_raeg_gtf)
  
  # ---------------------------------
  # Regiones genómicas
  # ---------------------------------
  coord_fau_espe = "s4:119661900-119665300"
  coord_actb_espe = "s13:2763300-2765500"
  
  coord_fau_raeg = "scaffold_m13_p_5:116490050-116493100"
  coord_actb_raeg = "scaffold_m13_p_7:71848400-71850750"
  
  # ---------------------------------
  # Generar los ocho SVG
  # ---------------------------------
  jobs = [
      ("01_eSpe_SRR6958995_FAU", "SRR6958995", coord_fau_espe, fau_gtfs["SRR6958995"]),
      ("02_eSpe_SRR6958995_ACTB", "SRR6958995", coord_actb_espe, actb_espe_gtf),
  
      ("03_eSpe_SRR6958999_FAU", "SRR6958999", coord_fau_espe, fau_gtfs["SRR6958999"]),
      ("04_eSpe_SRR6958999_ACTB", "SRR6958999", coord_actb_espe, actb_espe_gtf),
  
      ("05_eSpe_SRR6959000_FAU", "SRR6959000", coord_fau_espe, fau_gtfs["SRR6959000"]),
      ("06_eSpe_SRR6959000_ACTB", "SRR6959000", coord_actb_espe, actb_espe_gtf),
  
      ("07_rAeg_SRR18652228_FAU", "SRR18652228", coord_fau_raeg, fau_gtfs["SRR18652228"]),
      ("08_rAeg_SRR18652228_ACTB", "SRR18652228", coord_actb_raeg, actb_raeg_gtf),
  ]
  
  if not GGS.is_file():
      raise FileNotFoundError(f"No encontré ggsashimi: {GGS}")
  
  for name, run, coords, gtf in jobs:
      run_ggsashimi(
          bams[run],
          coords,
          gtf,
          SVG_DIR / name
      )
  
  print("\n¡Terminado!")
  print(f"GTF auxiliares: {GTF_DIR}")
  print(f"Gráficos SVG: {SVG_DIR}")
  PY
