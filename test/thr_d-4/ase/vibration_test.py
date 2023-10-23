from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.vibrations import Vibrations
from ase.vibrations import VibrationsData

# 水分子の定義
water = Atoms('H2O',
              positions=[(0.7680122628,0.5893398150,0.0000000000),
                         (-0.7680122628,0.5893398150,0.0000000000),
						 (0.0000000000,-0.0066796299,0.0000000000)],
              cell=[10,10,10])

input_data = {
        'control': { 'calculation' : 'scf',
                     'prefix' : 'water',
                     'pseudo_dir' : '/Users/genmasahiko/Program/pp',
                     'outdir' : 'data'
                    },
        'system' : { 'ecutwfc' : 40.0,
                     'ecutrho' : 320.0,
                     'ibrav' : 0,
                    },
        }

pseudopotentials = {'H' : 'H.pbe-rrkjus_psl.1.0.0.UPF',
                    'O' : 'O.pbe-n-rrkjus_psl.1.0.0.UPF'}

# Quantum ESPRESSO計算用の設定
calc = Espresso(
    kpts=(1, 1, 1),  # kポイントメッシュ
    pseudopotentials = pseudopotentials,
    tprnfor = True,
    input_data = input_data)

water.set_calculator(calc)

# 振動モード計算
vib = Vibrations(water, name='ase.vib', delta=0.01, nfree=2)
vib.run()


# ゼロ点振動エネルギーの計算
vib.summary()
