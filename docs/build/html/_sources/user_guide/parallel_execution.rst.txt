========
並列実行
========

並列実行について
================

大規模なネットワークの構築は、それだけで計算量と時間がかかる操作です。BMTKがそうであるように、bionet_liteも並列実行に対応しています。

基本的な使い方
==============

bionet_liteをインポートしたネットワーク構築コードを、mpirunなどのコマンドで実行するだけで自動的に並列実行が行われます。

.. code-block:: bash

   # 4プロセスで並列実行
   mpirun -n 4 python build_network.py

.. code-block:: bash

   # 16プロセスで並列実行
   mpirun -n 16 python build_network.py

コードの変更は不要
==================

並列実行のために特別なコードの変更は必要ありません。bionetと全く同じように、単にmpirunで実行するだけです：

.. code-block:: python

   # build_network.py
   from bionetlite import NeuliteBuilder as NetworkBuilder

   # 通常通りのコード
   net = NetworkBuilder('v1')
   net.add_nodes(...)
   net.add_edges(...)
   net.build()
   net.save_nodes(output_dir='network')
   net.save_edges(output_dir='network')

このファイルをmpirunで実行：

.. code-block:: bash

   mpirun -n 8 python build_network.py

並列実行のしくみ
================

bionet_liteは、bionetの並列実行機構をそのまま利用しています：

1. 各MPIプロセスがネットワークの一部を担当
2. ノードとエッジの生成が分散処理される
3. 最終的に結果が統合されてファイルに出力される

次のステップ
============

* :doc:`../advanced/allen_v1_model` - 大規模ネットワークの実例
