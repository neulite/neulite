====================
基本的な使い方
====================

bionet_liteを使用したネットワーク構築
=====================================

.. note::
   bionet_liteの使用にはBMTKのコード・データが必要です。多くの場合、既存のBMTKコードを再利用できます。

NeuliteBuilderの初期化
-----------------------

bionet_liteは、既存のBMTKコードでNeuliteBuilderをインポートするだけで利用可能です：

.. code-block:: python

   from bionetlite import NeuliteBuilder as NetworkBuilder

   # Create network
   net = NetworkBuilder(
       'network_name',
       convert_morphologies=True,  # Convert morphology files (default: True)
       convert_ion_channels=True   # Convert ion channel settings (default: True)
   )

パラメータの説明
^^^^^^^^^^^^^^^^

* ``network_name``: ネットワークの名前（文字列）
* ``convert_morphologies``: swcファイルをNeulite用に前処理するか（True/False）
* ``convert_ion_channels``: イオンチャネル設定ファイル（JSON）をNeulite用のCSVに変換するか（True/False）

.. note::
   ``convert_morphologies`` と ``convert_ion_channels`` はデフォルトでTrueに設定されています。2回目以降、Falseに設定することで処理をスキップできます。

ノードの追加
------------

既存のBMTKコードとノードの追加方法は同一です。

.. code-block:: python

   net.add_nodes(
       N=100,  # Number of neurons
       pop_name='Scnn1a',  # Population name
       model_type='biophysical',  # Model type
       model_template='ctdb:Biophys1.hoc',
       dynamics_params='472363762_fit.json',  # Ion channel parameters
       morphology='Scnn1a_473845048_m.swc',  # Morphology file
       ei='e'  # Excitatory ('e') or inhibitory ('i')
   )

サポートされているモデルタイプ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

現在、Bionet_liteならびにNeuliteは生物物理学的ニューロンのみに対応しています。Point neuronやVirtual cellは自動的に出力から除外されます。

* ✅ ``biophysical``: 形状を持つ生物物理学的ニューロン
* ❌ ``point_neuron``: ポイントニューロン（LIFなど）- 非対応
* ❌ ``virtual``: バーチャルニューロン - 非対応

エッジ（シナプス接続）の追加
----------------------------

基本的な接続の定義
^^^^^^^^^^^^^^^^^^

.. code-block:: python

   net.add_edges(
       source={'pop_name': 'Exc'},  # 送信元ポピュレーション
       target={'pop_name': 'Inh'},  # 受信先ポピュレーション
       connection_rule=5,  # 接続数または接続確率
       syn_weight=0.001,  # シナプス重み
       delay=2.0,  # 遅延（ms）
       dynamics_params='AMPA_ExcToInh.json',  # シナプスパラメータ
       model_template='exp2syn'  # シナプスモデル
   )

サポートされているシナプスモデル
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

現在、bionet_liteは以下のシナプスモデルに対応しています：

* ✅ ``exp2syn``: Double exponential synapse
* ❌ その他のシナプスモデル - 非対応

.. note::
   exp2synモデルの場合、configで設定されたディレクトリにあるシナプスモデルを自動的に適用します。

ネットワークの保存
------------------

.. code-block:: python

   # ネットワークのビルド
   net.build()

   # ノード情報の保存
   net.save_nodes(output_dir='network')

   # エッジ情報の保存
   net.save_edges(output_dir='network')

実行の流れ
==========

bionet_liteの実行フローは以下の通りです：

1. **NeuliteBuilderのインポート**: bionetの代わりにbionet_liteをインポート
2. **ネットワークの定義**: bionetと同じ方法でノードとエッジを定義
3. **ビルド**: ``net.build()`` でネットワークを構築
4. **保存**: ``save_nodes()`` と ``save_edges()`` でファイルを生成
5. **Neuliteでの実行**: 生成されたファイルをNeuliteカーネルで実行


制限事項
--------

bionet_liteは以下の機能に対応していません：

* Point neuron (LIF, etc.)
* Virtual neuron
* ファイルからのスパイク入力
* exp2syn以外のシナプスモデル

詳細は :doc:`../advanced/specification` を参照してください。

次のステップ
============

* :doc:`../api_reference/file_formats` - 生成されるファイルの詳細
* :doc:`configuration` - 設定ファイルの記述方法
* :doc:`parallel_execution` - 並列実行の方法
