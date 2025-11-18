==============
仕様と制限事項
==============

bionet_liteとNeuliteは、現在最小の実装と `富岳 <https://www.r-ccs.riken.jp/en/fugaku/>`_ 向けの報告が行われていますが、以下の仕様と制限事項があります。

サポートされている機能
======================

ニューロンモデル
----------------

* ✅ **Biophysical neurons**: 形状を持つ生物物理学的ニューロン（`Perisomatic model <https://brain-map.org/our-research/computational-modeling/perisomatic-biophysical-single-neuron-models>`_）
* ❌ **Point neurons**: ポイントニューロン（LIF、IAFなど） - 非対応
* ❌ **Virtual neurons**: バーチャルニューロン（入力用） - 非対応

シナプスモデル
--------------

* ✅ **exp2syn**: Double exponential synapse model
* ❌ その他のシナプスモデル - 非対応

ただし、exp2synであればconfigで設定されたディレクトリにあるシナプスモデルパラメータを利用可能です。

入力方式
--------

* ✅ **定常電流入力**: デフォルトでサポート
* ❌ **ファイルからのスパイク列入力**: 非対応（Neulite側での実装が必要）
* ❌ **Virtual cellによる入力**: 非対応

.. _形態モデル:

形態モデル
----------

現在bionet_liteは **Perisomatic model** に対応しています。

Perisomatic Modelとは
^^^^^^^^^^^^^^^^^^^^^

Perisomatic modelは、`Allen Cell Types Database <https://celltypes.brain-map.org/>`_ で採用されている生物物理学的ニューロンモデルの一種です。このモデルには以下の特徴があります：

**定義と構造:**

* **イオンチャネルの配置**: 主要なイオンチャネルは細胞体周辺に配置されます
* **詳細な形態**: 実際のニューロンから再構築された形態データを使用
* **簡略化された軸索**: 実際の軸索形態の代わりに、人工的な簡易軸索が付与されます

.. seealso::
   Allen Instituteの公式ドキュメント: `Perisomatic Single Neurons <https://brain-map.org/our-research/computational-modeling/perisomatic-biophysical-single-neuron-models>`_

制限事項の詳細
==============

1. ニューロンモデルの制限
-------------------------

bionetでは可能だが、bionet_liteで未対応のこと
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Point neuronとVirtual neuronは非対応**

BMTKやbionetには、形状を持たないpoint neuronや、実態は持たず入力元として利用されるVirtual neuronを定義することが可能ですが、bionet_liteならびにNeuliteはこれに対応していません。

bionet_liteは以下のニューロンタイプのみに対応：

* ✅ ``model_type='biophysical'``
* ❌ ``model_type='point_neuron'``
* ❌ ``model_type='virtual'``

2. 入力方式の制限
-----------------

**スパイク入力は非対応**

bionet_liteは入力方式のデフォルトとして、ニューロンへの定常電流入力のみを受け付けています。

BMTKはH5ファイルなどの形式で保存されたスパイク列をシミュレーションの入力として利用可能ですが、bionet_liteはこれに対応していません。これはNeuliteが定常電流入力のみが定義可能であることに由来します。

.. note::
   スパイク入力機能が必要な場合は、Neulite側での実装が必要です。

3. シナプスモデルの制限
-----------------------

**exp2syn以外は非対応**

シナプス接続については、現在は **double exponential synapse model (exp2syn)** にのみ対応しています。

exp2syn以外のシナプスモデルには非対応ですが、exp2synであればconfigで設定されたディレクトリにあるシナプスモデルを利用可能です。

4. シナプス接続位置の制限
-------------------------

**詳細な接続位置の手動指定は非対応**

現在のbionet_liteは、細胞間のシナプス接続位置を詳細に設定することはできません。

bionetはネットワーク構築時にシナプス接続のルール（特定のポピュレーションからポピュレーションへ、など）を設定することでネットワークを構築します。

しかし、実際にニューロン上のセクションaにおけるsegmentのbにおける距離xの位置、といった詳細度で接続位置が決定されるのは **bionet.simulatorモジュールが実行されたタイミング** となります。

現在のbionet_liteは、ネットワーク構築時に設定されたルールに該当する接続の位置を **ランダムに決定** します。

したがって：

* ❌ bionetと全く同じシナプス接続を取得することは現在できません
* ✅ 一度生成したNeulite用のbionet_liteによる出力は、本来シミュレーション実行時に決める詳細な位置情報を含み、確実に同じネットワークに対してシミュレーションを実行可能

5. データ型の制限
-----------------

**Delay値はint型のみ**

.. warning::
   BMTKのdelay値はdouble型ですが、**Neuliteはint型のみを受け付けます**。

   小数点以下のdelay値を使用する場合は、適切に変換する必要があります。

6. Perisomatic modelの注意点
----------------------------

bionet_liteが採用しているPerisomatic modelについては、:ref:`形態モデル` セクションを参照してください。

特に、実際の軸索形態が簡略化された人工軸索に置き換えられることに注意が必要です。

対応状況のまとめ
================

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - 機能
     - 対応状況
     - 備考
   * - Biophysical neurons
     - ✅ 対応
     - Perisomatic modelのみ
   * - Point neurons
     - ❌ 非対応
     -
   * - Virtual neurons
     - ❌ 非対応
     -
   * - 定常電流入力
     - ✅ 対応
     -
   * - スパイク列入力
     - ❌ 非対応
     - Neulite側実装が必要
   * - exp2synシナプス
     - ✅ 対応
     -
   * - その他のシナプス
     - ❌ 非対応
     -
   * - 並列実行
     - ✅ 対応
     -
   * - Allen Cell Types Databaseモデル
     - ✅ 対応
     -
   * - `SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ 形式
     - ✅ 対応
     - bionet経由

今後の発展
==========

bionet_lite、Neuliteは現在最小の実装と富岳向けの報告が行われていますが、今後の発展に伴いこれらの仕様は変更される可能性があります。

次のステップ
============

* :doc:`../tutorials/tutorial01_single_cell` - 実践的な使用例
