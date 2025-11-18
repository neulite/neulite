====================
よくある質問（FAQ）
====================

一般的な質問
============

bionet_liteとは何ですか？
--------------------------

bionet_liteは、BMTKのbionetモジュールを拡張し、新しい生物物理学的ニューロンシミュレータ「Neulite」のためのネットワークデータを生成するPythonパッケージです。

なぜbionet_liteが必要なのですか？
---------------------------------

bionet_liteは、Neuliteカーネルを使って生物物理学的ニューロンシミュレーションを実行するために必要です。

Neuliteは軽量で高性能なシミュレータで、スーパーコンピュータでの大規模シミュレーションに最適化されています。bionet_liteは、BMTKのネットワーク構築機能を活用しながら、Neulite用のネットワークファイルを生成します。

詳細は以下を参照してください：

* :doc:`neulite` - Neuliteカーネルについて
* :doc:`bionet_lite` - bionet_liteについて
* :doc:`architecture/design_and_implementation` - BMTKとBioNetの背景

bionetとbionet_liteの違いは何ですか？
--------------------------------------

主な違いは以下の通りです：

* **ターゲットシミュレータ**: bionet→NEURON、bionet_lite→Neulite
* **接続位置の決定**: bionetはシミュレータ実行時、bionet_liteはネットワーク構築時
* **出力ファイル**: bionet_liteはNeulite用のCSVファイルを追加生成

インストールと環境
==================

必要なPythonバージョンは？
--------------------------

Python 3.7以上が必要です。

詳細なセットアップ手順については、:doc:`getting_started/bionet_lite_setup` を参照してください。

どのOSで動作しますか？
----------------------

* Linux（推奨）
* macOS

スーパーコンピュータで実行できますか？
--------------------------------------

はい、Neuliteはスーパーコンピュータでの実行を想定して設計されています。特に `富岳 <https://www.r-ccs.riken.jp/en/fugaku/>`_ での実行実績があります。

一方で、bionet_liteはPythonの複数のライブラリに依存しているため、bionet_liteを用いたネットワーク構築はクラスタマシン等の環境で実行することを推奨しています。

使用方法
========

既存のbionetコードをbionet_liteに変換するには？
------------------------------------------------

import文を変更するだけです：

.. code-block:: python

   # Before change
   from bmtk.builder.networks import NetworkBuilder

   # After change
   from bionetlite import NeuliteBuilder as NetworkBuilder

残りのコードはそのまま使用できます。

Point neuronやVirtual neuronは使えますか？
-------------------------------------------

いいえ、現在bionet_liteはbiophysicalニューロンのみをサポートしています。

なぜexp2syn以外のシナプスモデルは使えないのですか？
---------------------------------------------------

現在のNeuliteの実装がexp2synのみをサポートしているためです。他のシナプスモデルのサポートは今後の開発課題です。

スパイク入力は使えますか？
--------------------------

現在、bionet_liteとNeuliteは定常電流入力のみをサポートしています。スパイク入力機能は、Neulite側での実装が必要です。

並列実行
========

並列実行の方法は？
------------------

mpirunを使用するだけです：

.. code-block:: bash

   mpirun -n 8 python build_network.py

コードの変更は不要です。詳細は :doc:`user_guide/parallel_execution` を参照してください。


トラブルシューティング
======================

"Module not found" エラーが出ます
----------------------------------

bionet_liteが正しくセットアップされているか確認してください：

.. code-block:: bash

   # Verify bionetlite.py exists
   ls bionetlite.py

セットアップされていない場合は、:doc:`getting_started/bionet_lite_setup` を参照してください。

接続位置がbionetと異なります
------------------------------

これは仕様です。bionet_liteは接続位置をネットワーク構築時にランダムに決定します。bionetとは異なるシード値を使用するため、完全に同一の接続位置にはなりません。

ただし、一度生成した<src>_<trg>_connection.csvを使用すれば、確実に同じネットワークでシミュレーションを実行できます。

delay値が整数しか受け付けられません
------------------------------------

これはNeuliteの現在の制限です。delay値はint型のみをサポートしています。小数点以下のdelay値が必要な場合は、時間刻み幅を調整するか、Neulite側の修正が必要です。



関連リソース
============

* :doc:`getting_started/index` - 入門ガイド
* :doc:`user_guide/basic_usage` - ユーザーガイド
* :doc:`advanced/specification` - 詳細な仕様
