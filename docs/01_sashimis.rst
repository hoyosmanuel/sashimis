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
