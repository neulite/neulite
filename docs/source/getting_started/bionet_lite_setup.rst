============================
bionet_liteのセットアップ
============================

このページでは、ネットワーク構築用のPythonパッケージである **bionet_lite** のセットアップ手順を説明します。

前提条件
========

以下のソフトウェアが必要です：

* Python 3.7以上

セットアップ手順
================

ステップ1: 仮想環境の作成（推奨）
----------------------------------

Python仮想環境の作成を推奨します：

.. code-block:: bash

   # Create virtual environment
   python -m venv neulite_env

   # Activate virtual environment (Linux/Mac)
   source neulite_env/bin/activate


ステップ2: BMTKとNEURONのインストール
--------------------------------------

bionet_liteは`Brain Modeling Toolkit (BMTK) <https://alleninstitute.github.io/bmtk/>`__と `NEURON <https://www.neuron.yale.edu/neuron/>`_ に依存しています。

**BMTKのインストール:**

.. code-block:: bash

   pip install bmtk


**NEURONのインストール:**

NEURON（生物物理学的ニューロンシミュレータ）をインストールします：

.. code-block:: bash

   pip install neuron

.. note::
   BMTK、NEURONについては公式サイトを参照してください：`BMTK <https://alleninstitute.github.io/bmtk/>`__ / `NEURON <https://www.neuron.yale.edu/neuron/>`__

ステップ3: 必要なパッケージのインストール
------------------------------------------

**scikit-learnのインストール:**

bionet_liteは軸索方向推定にPCA（主成分分析）を使用するため、scikit-learnが必要です：

.. code-block:: bash

   pip install scikit-learn

ステップ4: bionet_liteのダウンロード
-------------------------------------

bionet_liteのパッケージをダウンロードします：

.. code-block:: bash

   # Download file from download link
   wget [bionet_lite download URL]

   # Or use curl
   curl -L -O [bionet_lite download URL]


ステップ5: bionet_liteの配置
----------------------------

ダウンロードしたファイルを展開し、 ``bionetlite.py`` をネットワーク構築コードと同じディレクトリに配置します：

.. code-block:: bash

   # Extract tarball
   tar -xzf bionet_lite-X.Y.Z.tar.gz

   # Copy bionetlite.py to working directory
   cp bionet_lite-X.Y.Z/bionetlite.py /path/to/your/network/code/

**ディレクトリ構成：**

``bionetlite.py`` をネットワーク構築スクリプトと同じ階層に配置します：

.. code-block:: text

   your_project/
   ├── bionetlite.py           # Copy here
   ├── build_network.py        # Your network construction script
   └── config.json             # Configuration file (if needed)

この配置により、 ``build_network.py`` 内で以下のようにimportできます：

.. code-block:: python

   # build_network.py
   from bionetlite import NeuliteBuilder as NetworkBuilder

   net = NetworkBuilder('my_network')
   net.add_nodes(...)
   net.build()
   net.save_nodes(output_dir='network')

ステップ6: 動作確認
--------------------

配置が正しく完了したことを確認します：

.. code-block:: bash

   # Move to network construction code directory
   cd /path/to/your/network/code/

   # Verify that bionetlite.py exists
   ls bionetlite.py

.. code-block:: python

   # Start Python interpreter to verify
   python

   >>> from bionetlite import NeuliteBuilder
   >>> print("bionet_lite is ready to use!")

エラーが表示されなければ、セットアップは成功です。

基本的な使い方
==============

既存のbionetコードをbionet_liteに変換するのは非常に簡単です。import文を変更するだけで、残りのコードはそのまま使用できます。

importの変更
------------

**変更前（bionet）:**

.. code-block:: python

   from bmtk.builder.networks import NetworkBuilder

**変更後（bionet_lite）:**

.. code-block:: python

   from bionetlite import NeuliteBuilder as NetworkBuilder

ネットワークの構築
------------------

以降のコードはbionetと全く同じです：

.. code-block:: python

   # Create network
   net = NetworkBuilder('v1')

   # Add nodes
   net.add_nodes(
       N=1,
       pop_name='Scnn1a',
       model_type='biophysical',
       model_template='ctdb:Biophys1.hoc',
       dynamics_params='472363762_fit.json',
       morphology='Scnn1a_473845048_m.swc'
   )

   # Build network
   net.build()

   # Save to files
   net.save_nodes(output_dir='network')
   net.save_edges(output_dir='network')

出力ファイルの確認
------------------

bionet_liteは、``neulite/`` ディレクトリ内にNeulite用のファイルを生成します：

**ディレクトリ構造**

.. code-block:: text

   neulite/
   ├── v1_population.csv              # Population information
   ├── v1_v1_connection.csv           # Connection information
   └── data/
       ├── *.swc                       # Processed morphology files
       └── *.csv                       # Ion channel settings

.. note::
   ファイル名は ``NetworkBuilder`` で指定した名前（この例では ``'v1'``）に基づいて生成されます。

各ファイルの詳細：

**<network_name>_population.csv** - ニューロンのポピュレーション情報

.. code-block:: text

   #n_cell,n_comp,name,swc_file,ion_file

* ``#n_cell``: Number of cells in the population
* ``n_comp``: Number of compartments
* ``name``: Population name
* ``swc_file``: Path to morphology file
* ``ion_file``: Path to ion channel configuration file

**<src>_<trg>_connection.csv** - シナプス接続情報

.. code-block:: text

   #pre nid,post nid,post cid,weight,tau_decay,tau_rise,erev,delay,e/i

* ``#pre nid``: Presynaptic neuron ID
* ``post nid``: Postsynaptic neuron ID
* ``post cid``: Postsynaptic compartment ID
* ``weight``: Synaptic weight
* ``tau_decay``: Decay time constant (ms)
* ``tau_rise``: Rise time constant (ms)
* ``erev``: Reversal potential (mV)
* ``delay``: Transmission delay (ms)
* ``e/i``: Excitatory (e) or inhibitory (i)

オプション設定
==============

NeuliteBuilderは、以下のオプションパラメータを受け取ることができます：

.. code-block:: python

   net = NetworkBuilder(
       'v1',
       convert_morphologies=True,  # Whether to convert SWC files
       convert_ion_channels=True   # Whether to convert ion channel settings
   )

``convert_morphologies`` と ``convert_ion_channels`` はデフォルトでTrueです。2回目以降の実行では、これらをFalseに設定することで変換をスキップできます。

依存パッケージの詳細
====================

bionet_liteの主要な依存パッケージ
----------------------------------

bionet_liteは以下のPythonパッケージに依存しています：

.. code-block:: text

   bmtk>=1.0.0                    # Brain Modeling Toolkit
   neuron>=8.0.0                  # NEURON simulator (required)
   numpy>=1.19.0                  # Numerical computing
   pandas>=1.2.0                  # Data manipulation
   h5py>=3.0.0                    # HDF5 file support
   scikit-learn>=0.24.0           # PCA for axon direction estimation (required)
   mpi4py>=3.0.0                  # MPI parallel processing (optional)

**必須パッケージ:**

.. code-block:: bash

   pip install bmtk neuron scikit-learn

**オプションパッケージ:**

.. code-block:: bash

   # When using MPI parallel processing
   pip install mpi4py

次のステップ
============

* :doc:`neulite_build` - Neuliteカーネルのビルド

* :doc:`../tutorials/tutorial01_single_cell` - 実践的なチュートリアル
* :doc:`../user_guide/basic_usage` - 詳細な使い方

参考資料
========

* `NEURON公式サイト <https://www.neuron.yale.edu/neuron/>`__
* `BMTK Installation Guide <https://alleninstitute.github.io/bmtk/installation.html>`__
* `scikit-learn Documentation <https://scikit-learn.org/>`_
* `mpi4py Documentation <https://mpi4py.readthedocs.io/>`_
* Neulite公式サイト: https://numericalbrain.org/neulite/
