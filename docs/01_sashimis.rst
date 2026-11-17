I. Install
==================

.. code-block:: bash

  conda activate samtools_env
  python -m pip install --only-binary=:all: pysam==0.24.1
  cd ~/herramental/ggsashimi
  curl -fL -O https://raw.githubusercontent.com/guigolab/ggsashimi/master/ggsashimi.py

Me faltaba esta biblioteca de R ...meh!

.. code-block:: bash

  Rscript -e 'install.packages("svglite", repos="https://cloud.r-project.org")'  


II. Get the data
================

Rebuild the master Sashimi TSV directly from the original species-specific
StringTie2 GTF and TOGA BED12 files. Run the commands below in macOS Terminal,
one block at a time. The curation of transcript IDs is a manual input: this
manifest records the specific transcripts selected in the IGV figures.

Prerequisites
-------------

* The original files are under ``/Volumes/Expansion/IGV/<Species>/``.
* Each run has one matching ``<accession>*.gtf``.
* Each species has one matching ``*TOGA.bed``.
* Python 3.9+ is required; no extra Python packages are needed for the table.
* These steps do not modify original BAM, GTF, or BED files.

Step 1. Create the working directory
------------------------------------

.. code-block:: bash

   mkdir -p ~/Desktop/shashimis
   cd ~/Desktop/shashimis

Step 2. Record the manually curated transcript selection
---------------------------------------------------------

This input list contains 57 explicitly selected transcripts from 30 unique
species/run combinations. Each line is species, gene, run accession,
transcript ID, and the tissue label from the working table. Tissue labels
must be checked independently: the existing source list conflicts for
``Myotis_daubentonii / ERR13757853`` (Brain vs Heart).

Paste the entire block (including the terminating ``EOF``):

.. code-block:: bash

   cat > ~/Desktop/shashimis/selected_transcripts.tsv <<'EOF'
   Species	Gene	Run_accession	Transcript_ID	Tissue
   Eonycteris_spelaea	FAU	SRR6958995	STRG.3514.1	Mixed
   Eonycteris_spelaea	FAU	SRR6958999	STRG.4849.2	Mixed
   Eonycteris_spelaea	FAU	SRR6959000	STRG.4297.1	Mixed
   Eonycteris_spelaea	ZFP28	SRR6951024	STRG.5451.3	Mixed
   Myotis_daubentonii	KATNAL1	ERR13757851	STRG.3154.1	Brain
   Myotis_daubentonii	NCL	ERR13757850	STRG.11413.1	Kidney
   Myotis_daubentonii	NCL	ERR13757851	STRG.10405.1	Brain
   Myotis_daubentonii	NCL	ERR13757851	STRG.10405.2	Brain
   Myotis_daubentonii	NCL	ERR13757851	STRG.10405.3	Brain
   Myotis_daubentonii	NCL	ERR13757852	STRG.11816.1	Lung
   Myotis_daubentonii	TMPO	ERR13757850	STRG.22537.1	Kidney
   Myotis_daubentonii	TMPO	ERR13757851	STRG.20121.1	Brain
   Myotis_daubentonii	TMPO	ERR13757852	STRG.23070.1	Lung
   Myotis_daubentonii	TMPO	ERR13757853	STRG.15021.1	Brain
   Myotis_daubentonii	ZNF24	ERR13757853	STRG.8426.2	Heart
   Myotis_myotis	NCL	SRR18652219	STRG.3566.1	Brain
   Myotis_myotis	TMPO	SRR18652219	STRG.7092.1	Brain
   Myotis_mystacinus	NCL	ERR13757860	STRG.9638.1	
   Myotis_mystacinus	NCL	ERR13757861	STRG.11330.1	Lung
   Myotis_mystacinus	NCL	ERR13757861	STRG.11330.2	Lung
   Myotis_mystacinus	NCL	ERR13757862	STRG.13487.1	Spleen
   Myotis_mystacinus	NCL	ERR13757862	STRG.13487.2	Spleen
   Myotis_mystacinus	NCL	ERR13757863	STRG.14317.1	Mixed
   Myotis_mystacinus	NCL	ERR13757863	STRG.14317.2	Mixed
   Myotis_mystacinus	NCL	ERR13757864	STRG.13567.2	Kidney
   Myotis_mystacinus	NCL	ERR13757864	STRG.13567.3	Kidney
   Myotis_mystacinus	TMPO	ERR13757861	STRG.21823.1	Lung
   Myotis_mystacinus	TMPO	ERR13757863	STRG.26900.1	Mixed
   Plecotus_auritus	AP1AR	ERR13757871	STRG.8744.3	Mixed
   Plecotus_auritus	AP1AR	ERR13757872	STRG.11209.3	Mixed
   Plecotus_auritus	NCL	ERR13757870	STRG.6157.1	Mixed
   Plecotus_auritus	NCL	ERR13757870	STRG.6157.2	Mixed
   Plecotus_auritus	NCL	ERR13757870	STRG.6157.3	Mixed
   Plecotus_auritus	NCL	ERR13757871	STRG.6561.1	Mixed
   Plecotus_auritus	NCL	ERR13757872	STRG.8419.1	Mixed
   Plecotus_auritus	NCL	ERR13757873	STRG.8375.3	Mixed
   Plecotus_auritus	NCL	ERR13757873	STRG.8375.4	Mixed
   Plecotus_auritus	NCL	ERR13757873	STRG.8375.5	Mixed
   Rhinolophus_ferrumequinum	FBXO3	SRR18652180	STRG.3111.1	Brain
   Rhinolophus_ferrumequinum	ZKSCAN7	SRR18652180	STRG.4351.1	Brain
   Rhinolophus_hipposideros	FBXO3	ERR13757944	STRG.9249.2	Lung
   Rhinolophus_hipposideros	FBXO3	ERR13757945	STRG.6028.1	Brain
   Rhinolophus_hipposideros	FBXO3	ERR13757946	STRG.11571.3	Kidney
   Rhinolophus_hipposideros	ZKSCAN7	ERR13757944	STRG.13078.1	Lung
   Rhinolophus_hipposideros	ZKSCAN7	ERR13757945	STRG.8409.1	Brain
   Rhinolophus_hipposideros	ZKSCAN7	ERR13757946	STRG.16281.1	Kidney
   Rhinolophus_hipposideros	ZKSCAN7	ERR13757947	STRG.12323.1	Heart
   Rhinolophus_sinicus	FBXO3	SRR17418361	STRG.6148.1	Mixed
   Rhinolophus_sinicus	ZKSCAN7	SRR17418361	STRG.9548.1	Mixed
   Rousettus_aegyptiacus	FAU	SRR18652228	STRG.3989.1	Testes
   Rousettus_aegyptiacus	FAU	SRR18652228	STRG.3989.3	Testes
   Rousettus_aegyptiacus	ZFP28	SRR18652228	STRG.7996.1	Testes
   Vespertilio_murinus	NCL	ERR13757845	STRG.9901.1	Stomach
   Vespertilio_murinus	NCL	ERR13757846	STRG.13689.1	Brain
   Vespertilio_murinus	NCL	ERR13757847	STRG.15933.1	Kidney
   Vespertilio_murinus	NCL	ERR13757848	STRG.16627.1	Lung
   Vespertilio_murinus	NCL	ERR13757849	STRG.11029.1	Liver
   EOF


Step 3. Create the extractor script
------------------------------------

This script reads source files, matches exact transcript IDs, extracts
transcript spans, exon coordinates and TPM, locates ACTB from TOGA BED12,
collapses identical ACTB structural duplicates, converts GTF to 0-based
half-open coordinates, and emits one ACTB control per selected experiment.
It fails if an input is missing or the ACTB annotation is ambiguous.

.. code-block:: bash

   cat > ~/Desktop/shashimis/build_sashimi_master.py <<'PY'
   #!/usr/bin/env python3
   """Rebuild the Sashimi master TSV from a curated selection and source GTF/TOGA BED12.
   
   Standard library only. No BAMs or source annotations are modified.
   """
   import csv
   import re
   import sys
   from pathlib import Path
   from collections import defaultdict
   
   ROOT = Path('/Volumes/Expansion/IGV')
   WORK = Path.home() / 'Desktop' / 'shashimis'
   SELECTION = WORK / 'selected_transcripts.tsv'
   OUTPUT = WORK / 'Base_maestra_ggsashimi_reconstruida.tsv'
   
   FIELDS = [
       'Species','Gene','Role','Chimera_type','Run_accession','Tissue',
       'Transcript_ID','Assembly','Scaffold','Start_0based','End_exclusive',
       'Strand','Annotation_source','Source_file','TPM','Status','Notes'
   ]
   PARALOG = {'FAU','FBXO3','NCL','TMPO'}
   ISOFORM = {'AP1AR','KATNAL1','ZFP28','ZKSCAN7','ZNF24'}
   
   
   def choose_one(folder, pattern):
       found = sorted(p for p in folder.glob(pattern) if p.is_file() and not p.name.startswith('._'))
       if len(found) != 1:
           raise ValueError(f'{folder}: expected exactly one {pattern}, got {len(found)}: {[p.name for p in found]}')
       return found[0]
   
   
   def read_selection():
       with SELECTION.open(newline='') as fh:
           reader = csv.DictReader(fh, delimiter='\t')
           required = {'Species','Gene','Run_accession','Transcript_ID','Tissue'}
           if not required.issubset(reader.fieldnames or []):
               raise ValueError(f'Selection must contain columns: {sorted(required)}')
           records = list(reader)
       keys = [(r['Species'],r['Gene'],r['Run_accession'],r['Transcript_ID']) for r in records]
       if len(keys) != len(set(keys)):
           raise ValueError('Duplicated species/gene/run/transcript in selection')
       for rec in records:
           if not all(rec[k].strip() for k in ('Species','Gene','Run_accession','Transcript_ID')):
               raise ValueError(f'Incomplete selection row: {rec}')
           if rec['Gene'] not in PARALOG | ISOFORM:
               raise ValueError(f'Unknown gene: {rec["Gene"]}')
       return records
   
   
   def selected_transcript(gtf, transcript_id):
       needle = re.compile(r'transcript_id\s+"' + re.escape(transcript_id) + r'"\s*;')
       txs, exons = [], []
       with gtf.open() as fh:
           for raw in fh:
               if raw.startswith('#') or not needle.search(raw):
                   continue
               f = raw.rstrip('\n').split('\t')
               if len(f) < 9:
                   raise ValueError(f'Malformed GTF in {gtf}')
               if f[2] == 'transcript':
                   txs.append(f)
               elif f[2] == 'exon':
                   exons.append(f)
       if len(txs) != 1 or not exons:
           raise ValueError(f'Expected 1 transcript and >=1 exon for {transcript_id} in {gtf}: found {len(txs)} and {len(exons)}')
       t = txs[0]
       if any(e[0] != t[0] or e[6] != t[6] or int(e[3]) < int(t[3]) or int(e[4]) > int(t[4]) for e in exons):
           raise ValueError(f'Exons inconsistent with {transcript_id} in {gtf}')
       attr = re.search(r'\bTPM\s+"([^"]+)"', t[8])
       if not attr:
           raise ValueError(f'Missing TPM for {transcript_id} in {gtf}')
       exons.sort(key=lambda e: (int(e[3]), int(e[4])))
       exon_note = '; '.join(f'{e[3]}-{e[4]}' for e in exons)
       return t[0], int(t[3])-1, int(t[4]), t[6], attr.group(1), exon_note
   
   
   def actb_from_bed(bed):
       structs = {}
       with bed.open() as fh:
           for raw in fh:
               f = raw.rstrip('\n').split('\t')
               if len(f) < 12 or '#ACTB#' not in f[3]:
                   continue
               key = tuple(f[:3] + f[5:12])
               structs[key] = f
       if len(structs) != 1:
           raise ValueError(f'{bed}: expected one unique ACTB BED12 structure, found {len(structs)}')
       f = next(iter(structs.values()))
       start, end = int(f[1]), int(f[2])
       count = int(f[9])
       sizes = [int(x) for x in f[10].rstrip(',').split(',')]
       offsets = [int(x) for x in f[11].rstrip(',').split(',')]
       if not (len(sizes) == len(offsets) == count) or start >= end:
           raise ValueError(f'Invalid ACTB BED12 structure in {bed}')
       for offset, size in zip(offsets, sizes):
           if offset < 0 or size <= 0 or start + offset + size > end:
               raise ValueError(f'Invalid ACTB exon boundary in {bed}')
       return f[0], start, end, f[5], count
   
   
   def main():
       selections = read_selection()
       if not ROOT.is_dir():
           raise FileNotFoundError(ROOT)
       records, controls = [], {}
       tissue_by_run = defaultdict(set)
       for r in selections:
           species, gene, run, tid = (r[k].strip() for k in ('Species','Gene','Run_accession','Transcript_ID'))
           folder = ROOT / species
           gtf = choose_one(folder, f'{run}*.gtf')
           chrom, start, end, strand, tpm, exons = selected_transcript(gtf, tid)
           bed = choose_one(folder, '*TOGA.bed')
           assembly = bed.name.removesuffix('.TOGA.bed')
           target = dict.fromkeys(FIELDS, '')
           target.update(Species=species, Gene=gene, Role='Target',
               Chimera_type='Paralog' if gene in PARALOG else 'Isoform',
               Run_accession=run, Tissue=r['Tissue'].strip(), Transcript_ID=tid,
               Assembly=assembly, Scaffold=chrom, Start_0based=str(start),
               End_exclusive=str(end), Strand=strand,
               Annotation_source='StringTie2 GTF', Source_file=gtf.name, TPM=tpm,
               Status='GTF verified; BAM pending',
               Notes=f'{len(exons.split("; "))} exons; exon positions 1-based inclusive: {exons}')
           records.append(target)
           tissue_by_run[(species, run)].add(r['Tissue'].strip())
           control_key = species, run
           if control_key not in controls:
               actb_chrom, actb_start, actb_end, actb_strand, actb_exon_count = actb_from_bed(bed)
               control = dict.fromkeys(FIELDS, '')
               control.update(Species=species, Gene='ACTB', Role='Control',
                   Run_accession=run, Tissue=r['Tissue'].strip(),
                   Assembly=assembly, Scaffold=actb_chrom,
                   Start_0based=str(actb_start), End_exclusive=str(actb_end),
                   Strand=actb_strand, Annotation_source='TOGA BED12',
                   Source_file=bed.name, Status='ACTB BED12 verified; BAM pending',
                   Notes=f'{actb_exon_count} exons; identical duplicate BED12 records collapsed')
               controls[control_key] = control
       for key, tissues in tissue_by_run.items():
           nonempty = {t for t in tissues if t}
           if len(nonempty) > 1:
               controls[key]['Tissue'] = ''
               controls[key]['Status'] = 'ACTB BED12 verified; tissue conflict; BAM pending'
               controls[key]['Notes'] += f'; REVIEW: conflicting tissues in selection: {sorted(nonempty)}'
               print(f'WARNING: {key[0]} {key[1]} has conflicting tissue labels: {sorted(nonempty)}')
           elif len(nonempty) == 1:
               controls[key]['Tissue'] = next(iter(nonempty))
       combined = records + list(controls.values())
       combined.sort(key=lambda r: (r['Species'], r['Gene'], r['Run_accession'], r['Transcript_ID']))
       with OUTPUT.open('w', newline='') as fh:
           writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
           writer.writeheader()
           writer.writerows(combined)
       print(f'OK: {len(records)} target transcripts, {len(controls)} ACTB controls, {len(combined)} total rows')
       print(f'Output: {OUTPUT}')
   
   
   if __name__ == '__main__':
       try:
           main()
       except (OSError, ValueError) as err:
           sys.exit(f'ERROR: {err}')
   PY

Step 4. Run the extraction
--------------------------

.. code-block:: bash

   conda activate samtools_env
   python ~/Desktop/shashimis/build_sashimi_master.py

Expected output: ``57 target transcripts, 30 ACTB controls, 87 total rows``.
Any conflicting tissue annotation is reported as a warning and the ACTB
tissue field is left blank until checked.

Step 5. Inspect and validate the resulting TSV
-----------------------------------------------

.. code-block:: bash

   head -n 4 ~/Desktop/shashimis/Base_maestra_ggsashimi_reconstruida.tsv
   wc -l ~/Desktop/shashimis/Base_maestra_ggsashimi_reconstruida.tsv

The TSV has 88 lines including the header, and 17 columns. Verify counts:

.. code-block:: bash

   python - <<'PY'
   import csv
   from pathlib import Path
   p = Path.home() / 'Desktop/shashimis/Base_maestra_ggsashimi_reconstruida.tsv'
   with p.open(newline='') as f:
       rows = list(csv.DictReader(f, delimiter='\t'))
   targets = [r for r in rows if r['Role'] == 'Target']
   controls = [r for r in rows if r['Role'] == 'Control']
   assert len(targets) == 57, len(targets)
   assert len(controls) == 30, len(controls)
   assert len(rows) == 87, len(rows)
   assert len({(r['Species'],r['Run_accession']) for r in controls}) == 30
   assert {(r['Species'],r['Run_accession']) for r in targets} == {(r['Species'],r['Run_accession']) for r in controls}
   print('PASS: 57 targets, 30 ACTB controls, 87 records')
   PY

Cautions
--------

* This reconstructs the **selected, verified** records. The historical
  ``v024`` spreadsheet also included unselected pending placeholder records
  and orphan ACTB species-level records. Therefore this regenerated table
  has 87 rather than 94 rows; those additional rows are intentionally excluded.
* Source GTF provides transcript annotations, not direct proof that splice
  junctions are supported by independent molecules.
* The BAM headers and local junction support must be validated separately.
* Transcript IDs such as STRG.4297.1 are scoped to each individual GTF run.
