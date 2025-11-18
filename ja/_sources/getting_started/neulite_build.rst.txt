==========================
Neuliteカーネルのビルド
==========================

このページでは、シミュレーション実行用のCプログラムである **Neuliteカーネル** のビルド手順を説明します。

前提条件
========

以下のソフトウェアが必要です：

* C17対応のCコンパイラ（gcc 7.0以上推奨）
* make

ビルド手順
==========

ステップ1: ソースコードのダウンロード
--------------------------------------

Neuliteカーネルのソースコードをダウンロードします：

.. code-block:: bash

   # Download file from download link
   wget [Neulite kernel download URL]

   # Or use curl
   curl -L -O [Neulite kernel download URL]


ステップ2: ソースコードの展開
------------------------------

ダウンロードしたファイルを展開します：

.. code-block:: bash

   # Extract tarball
   tar -xzf neulite_kernel-X.Y.Z.tar.gz
   cd neulite_kernel-X.Y.Z

ステップ3: ビルド
------------------

Makefileを使用してビルドします：

.. code-block:: bash

   # Build
   make

これにより、実行ファイル ``nl`` が生成されます。

ステップ4: 実行
---------------

ビルドされた ``nl`` を実行します。引数としてBionet_liteが生成したポピュレーションファイル（例：``<network_name>_population.csv`` ）シナプス接続ファイル（ 例：``<src>_<trg>_connection.csv`` ）を指定します：

.. code-block:: bash

   # Run Neulite kernel
   ./nl <network_name>_population.csv <src>_<trg>_connection.csv

MPI並列実行する場合：

.. code-block:: bash

   # MPI parallel execution
   mpirun -np 4 ./nl <network_name>_population.csv <src>_<trg>_connection.csv

ビルドの確認
============

Neuliteカーネルが正しくビルドされたことを確認します：

.. code-block:: bash

   # Verify executable exists
   ls -l nl

   # Confirm file is executable
   file nl

トラブルシューティング
======================

コンパイルエラー
----------------

**問題**: ``implicit declaration of function``

**解決策**: C17対応のコンパイラを使用してください

.. code-block:: bash

   # Check gcc version
   gcc --version

   # gcc 7.0+ required. Update if older
   sudo apt-get install gcc-9



実行ファイルが見つからない
--------------------------

**問題**: ``./nl: No such file or directory``

**解決策**: ビルドディレクトリで実行していることを確認してください

.. code-block:: bash

   # Move to build directory
   cd /path/to/neulite_kernel-X.Y.Z

   # Verify file exists
   ls -l nl

次のステップ
============

Neuliteカーネルのビルドが完了したら、以下のチュートリアルに進んでください：

* :doc:`../tutorials/tutorial01_single_cell` - 実践的なチュートリアル
* :doc:`../user_guide/basic_usage` - 基本的な使い方
