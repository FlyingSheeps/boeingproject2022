import re
import control.matlab as mat
from matplotlib import pyplot as plt
import numpy as np
plt.rc('font', size=12)

def parse_text(text):
    parsed_data = []
    pattern = r"Calculation for control position ([+-]?\d+\.\d+)"
    matches = re.findall(pattern, text)
    
    for match in matches:
        control_position = float(match)
        control_data = {}
        
        pattern = rf"Calculation for control position {re.escape(match)}([\s\S]+?)_____Finished operating point calculation for control position {re.escape(match)}"
        section_match = re.search(pattern, text)
        if section_match:
            section_text = section_match.group(1)
            
            # Parse control derivatives
            pattern = r"Control derivatives([\s\S]+?)\n\n"
            match = re.search(pattern, section_text)
            if match:
                control_derivatives_text = match.group(1)
                control_derivatives = {}
                sub_pattern = r"(\w+)=(\s*[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
                sub_matches = re.findall(sub_pattern, control_derivatives_text)
                for sub_match in sub_matches:
                    key = sub_match[0]
                    value = float(sub_match[1])
                    control_derivatives[key] = value
                control_data['control_derivatives'] = control_derivatives
            
            # Parse state matrices
            pattern = r"Longitudinal state matrix([\s\S]+?)\n\n"
            match = re.search(pattern, section_text)
            if match:
                state_matrix_text = match.group(1)
                state_matrix = []
                rows = state_matrix_text.strip().split('\n')
                for row in rows:
                    if row == "       Lateral state matrix":
                        break
                    elements = [float(e) for e in row.split()]
                    state_matrix.append(elements)
                control_data['state_matrices'] = {'Longitudinal': state_matrix}
            
            pattern = r"Lateral state matrix([\s\S]+?)\n\n"
            match = re.search(pattern, section_text)
            if match:
                state_matrix_text = match.group(1)
                state_matrix = []
                rows = state_matrix_text.strip().split('\n')
                for row in rows:
                    elements = [float(e) for e in row.split()]
                    state_matrix.append(elements)
                control_data['state_matrices']['Lateral'] = state_matrix
            
            # Parse control matrices
            pattern = r"Longitudinal control matrix([\s\S]+?)\n\n"
            match = re.search(pattern, section_text)
            if match:
                control_matrix_text = match.group(1)
                control_matrix = []
                rows = control_matrix_text.strip().split('\n')
                for row in rows:
                    elements = [float(e) for e in row.split()]
                    control_matrix.append(elements)
                control_data['control_matrices'] = {'Longitudinal': control_matrix}
            
            pattern = r"Lateral control matrix([\s\S]+?)\n\n"
            match = re.search(pattern, section_text)
            if match:
                control_matrix_text = match.group(1)
                control_matrix = []
                rows = control_matrix_text.strip().split('\n')
                for row in rows:
                    elements = [float(e) for e in row.split()]
                    control_matrix.append(elements)
                control_data['control_matrices']['Lateral'] = control_matrix
            
            parsed_data.append(control_data)
    
    return parsed_data

# テキストデータを読み込むなどの前処理が必要な場合はここで行ってください
with open("xflr5.txt", "r") as file:
    text_data = file.read()


# テキストデータからデータをパースする
parsed_data = parse_text(text_data)

# 各データの取得例を表示する
for data in parsed_data:
    print("Control Derivatives:", data['control_derivatives'])
    print("Longitudinal State Matrix:", data['state_matrices']['Longitudinal'])
    print("Lateral State Matrix:", data['state_matrices']['Lateral'])
    print("Longitudinal Control Matrix:", data['control_matrices']['Longitudinal'])
    print("Lateral Control Matrix:", data['control_matrices']['Lateral'])
    print()

#状態空間表現から伝達関数を計算して，全部表示する
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for data in parsed_data:
    A = data['state_matrices']['Longitudinal']
    B = data['control_matrices']['Longitudinal']
    C = [0,0,0,1]

    sys = mat.ss(A,B,C,0)
    sys = mat.tf(sys)

    Kd = 7/100
    Ki = 0.8
    K = 0.6

    #サーボモデル
    omega_a = 2*3.14*3
    zeta_a = 2*0.8*omega_a
    servo = mat.tf([omega_a**2],[1,zeta_a,omega_a**2])
    #制御モデル
    controltf = mat.tf([Kd,0],[1]) + mat.tf([1],[1]) + mat.tf([Ki],[1,0])
    #システム全体
    sysall = sys*servo*controltf

    mag, phase, omega = mat.bode(sysall, omega_limits=[0.1,100],plot=False)
    # 位相に180度を加えます
    phase = phase*180/np.pi + 180
    mag = 20*np.log10(mag)
    
    # ゲインプロット
    ax1.semilogx(omega, mag)
    
    # 位相プロット
    ax2.semilogx(omega, phase)

# 軸ラベルを設定します
ax2.set_xlabel('Frequency [rad/s]')
ax1.set_ylabel('Gain [dB]')
ax2.set_ylabel('Phase [deg]')
ax2.set_yticks([30,0,-30,-60,-90,-120,-150,-180,-210,-240,-270])
ax1.grid(which="both")
ax2.grid(which="both")
plt.savefig("bodes.png")
plt.savefig("robust_W.eps")
plt.clf()

#根軌跡の計算
