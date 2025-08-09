# WASA2026
## Contents
全部で5種類あります。
1. WingProfileManager
> リブデータ図面自動生成プログラム
> 使用方法はフォルダ内のpdfを参照すること
2. Windmize
> WASAのOBである森田直人さんが製作したプログラムをMATLAB2024bで動くように修正したもの
3. TR-797
> 非平面翼の最適設計-揚力と翼根曲げモーメントを与えた時の最小誘導抵抗-の数値計算部分をMATLABで実装
> 論文中では非平面翼の計算を行なっていますが，ここに置いているスクリプトは平面翼です． スクリプト中のGeometric Conditionを変更することで非平面翼に対応
> 設計には直接使えないので参考程度に
4. Larrabee
> Larrabeeの方法によるプロペラ設計プログラム
5. Divergence
> ダイバージェンス速度解析プログラム
### Larrabeeフォルダ内容
#### Larrabee_main.m
- 実行すべきmファイル。ここから以下の各々のサブルーチンを実行する。
#### Larrabee_input.m
- プロペラの設定値を決める。ここの値を変えることで様々な条件のプロペラが設計できる。
#### readXFLR.m
- XFLR５で解析したairfoilフォルダ内のtxtフィイルを読み込み、 data_matという多次元配列に代入する。 例えば、'pelafoil_T1_Re0.01_M0.00_N9.0.txt'のようなファイルが必要。 それぞれの翼型に変えたければ中身のfoil_nameの値を変える。 レイノルズ数の数や値を変えたければRelistを変える。レイノルズ数は使用する領域より幅広く取ってやるとエラーが少なくなる。
#### Larrabee_airfoil.m
- Larrabee_airfoil_ini.mではCd（抗力係数）とalpha（迎え角）を自分で決めるが、readXFLR.mで翼型の解析データが読み込めればCdとalphaをCl（揚力係数）とRe（レイノルズ数）の2変数関数と置いて、(Cd = f(Re,Cl), alpha = f(Re,Cl))自分で決めたClの値と前回の計算から持ってきたReからCd、alphaを求める。
#### Larrabee_airfoil_ini.m
- Reの範囲がreadXFLRで読み込んだ値を越える場合は外挿ができないためにエラーが出るので、こちらを用いる。自分でCdとalphaの値を指定する。
#### Larrabee_culc.m
- 計算部分
#### Larrabee_output_fig.m
- 図として出力する。出力された図はresultフォルダに保存される。
#### Larrabee_output_result.m
- 計算された結果をresultフォルダのresult.txtに保存する。
#### Larrabee_output_cad.m
- 設計されたプロペラを３DCAD上で表示するためにCADデータをcsvなどで保存する。 SolidWorksに対応予定。翼型はairfoil内にある翼型.datファイルを読み込む。 未実装。
#### Re_lookup.m
- Larrabee_airfoil.mでCd、 alphaの値を補完するための関数。 interp1,interp2関数を用いて補完している。要改良。
### Divergenceフォルダ内容
#### wing.csvの構成
wing.csvはヘッダー行とそれに対応する値が列として入っているpandasで使用する前提のcsvファイルです． wing.csvを書き換える場合は以下の手順で行ってください．
1. まずspanにmm単位の翼根から翼端までの代表点位置を記述する．
2. その代表点位置と同じ行に代表点位置におけるねじり剛性などを入力する．
3. 空力係数を計算した際の釣り合い状態の飛行速度をU0に入力する．
4. 不等間隔のデータでも補完して等間隔に修正されて使用されます．これはscipy.interpolate.interp1dの許容するデータであれば問題ないです．

入力するデータの種類と単位は以下の通りです． MACは%ではなく割合で記入してください． U0は最初のセルのみ入力されていれば十分です．

- span：スパン方向位置 (mm)
- mass：区間の重量 (kg)
- EI：上下曲げ剛性 (Pa m^4)
- GIp：ねじり剛性 (Pa m^4)
- c：コード長 (mm)
- T.C.：ねじり中心or桁位置 (MAC)
- Cm：モーメント係数 (-)
- CL：揚力係数 (-)
- U0：入力した空力係数を計算した飛行速度 (m/s)
