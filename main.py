import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

massa_eletron = 9.11E-31
massa_proton = 1.67E-27
cte_plank_ev = 4.136E-15
cte_plank_joule = 6.626E-34
pi = math.pi
velocidade_luz = 3E8

def imprime_tela_inicial():
    print("\n***********************************************************************************")
    print("**************** SIMULADOR CONFINAMENTO DE PARTICULAS EM UMA CAIXA ****************")
    print("***********************************************************************************")

    print("\nAutoria: Mariane S. Carvalho\tTurma: 610")

    print("\n--------------------------------------------------------------------------------------------------------------\n")
    print("Este projeto tem como objetivo: \n\n\t- Simular o confinamento de particulas quanticas levando a quantizacao dos niveis de energia da mesma; \n\t- Entender os conceitos de salto quantico, variacao de niveis de energia atraves da emissao e absorção de fotons; \n\t- Estudar a função de onda quantica independente do tempo e a funcao de distribuicao de probabilidade |Ψ(x,t)|²")
    print()

def menu_nav():
    print("\n--------------------------------------------------------------------------------------------------------------")
    print("\nSelecione a entrada:\n\n\t1 - Determinar funcao de onda quantica e demais parametros;\n\t2 - Determinar parametros da particula e caixa a partir da funcao de onda\n\t3 - Encerrar programa")

def calcula_parametros_funcao_onda(largura_caixa, num_quantico):
    amplitude = math.sqrt((2) / (largura_caixa))
    num_onda = (num_quantico * pi) / (largura_caixa)

    return amplitude, num_onda

def plota_graficos_func_onda(amplitude_ni, amplitude_nf, num_onda_ni, num_onda_nf, largura_caixa, n_inicial, n_final):
    x = np.linspace(0, (largura_caixa), 150)

    fig, eixo = plt.subplots(1,2, sharey = True, figsize=(12,6))

    eixo[0].set_xlabel("x (m)", fontsize=16)
    eixo[1].set_xlabel("x (m)", fontsize=16)
    eixo[0].set_ylabel("Ψ₁", fontsize=16)
    eixo[1].set_ylabel("Ψ₂", fontsize=16)

    funcao_onda_ni = amplitude_ni * np.sin(num_onda_ni * x)
    funcao_onda_nf = amplitude_nf * np.sin(num_onda_nf * x)

    plot_func_ni, = eixo[0].plot(x, funcao_onda_ni)
    plot_func_nf, = eixo[1].plot(x, funcao_onda_nf)

    eixo[0].set_title(f"n = {n_inicial}")
    eixo[1].set_title(f"n = {n_final}")

    plt.savefig("graficos/grafico_funcao_onda.pdf")

    plt.close()

def plota_graficos_distr_probabilidade(amplitude_ni, amplitude_nf, num_onda_ni, num_onda_nf, largura_caixa, n_inicial, n_final):
    x = np.linspace(0, (largura_caixa), 150)

    fig, eixo = plt.subplots(1,2, sharey = True, figsize=(12,6))

    eixo[0].set_xlabel("x (m)", fontsize=16)
    eixo[1].set_xlabel("x (m)", fontsize=16)
    eixo[0].set_ylabel("|Ψ₁|²", fontsize=16)
    eixo[1].set_ylabel("|Ψ₂|²", fontsize=16)

    func_distr_probab_ni = (amplitude_ni ** 2) * (np.sin(num_onda_ni * x) ** 2)
    func_distr_probab_nf = (amplitude_nf ** 2) * (np.sin(num_onda_nf * x) ** 2)

    plot_func_ni, = eixo[0].plot(x, func_distr_probab_ni)
    plot_func_nf, = eixo[1].plot(x, func_distr_probab_nf)

    eixo[0].set_title(f"n = {n_inicial}")
    eixo[1].set_title(f"n = {n_final}")

    plt.savefig("graficos/grafico_distr_probab.pdf")

    plt.close()
    
def calcula_energia_particula(largura_caixa, num_quantico, tipo_particula):
    energia_joule = (math.pow(num_quantico, 2) * math.pow(cte_plank_joule, 2)) / (8 * math.pow(largura_caixa, 2) * massa_eletron)

    # se particula confinada for proton
    if tipo_particula == 2:
        energia_joule = (math.pow(num_quantico, 2) * math.pow(cte_plank_joule, 2)) / (8 * math.pow(largura_caixa, 2) * massa_proton)

    return energia_joule

def calcula_energia_foton(num_quantico_inicial, num_quantico_final, tipo_particula, largura_caixa):
    # se n_inicial > n_final: foton emitido
    # se n_inicial < n_final: foton absorvido
    energia_inicial = calcula_energia_particula(largura_caixa, num_quantico_inicial, 1)
    energia_final = calcula_energia_particula(largura_caixa, num_quantico_final, 1)

    if tipo_particula == 2:
        energia_inicial = calcula_energia_particula(largura_caixa, num_quantico_inicial, 2)
        energia_final = calcula_energia_particula(largura_caixa, num_quantico_final, 2)

    energia_foton = abs(energia_final - energia_inicial)

    return energia_foton

def calcula_comprimento_onda_freq_foton(energia_foton):
    comprimento_onda = (cte_plank_joule * velocidade_luz) / energia_foton

    freq = energia_foton / cte_plank_joule

    return comprimento_onda, freq

def calcula_velocidade_particula(energia_particula, massa_particula):
    velocidade = math.sqrt((2 * energia_particula) / massa_particula)

    return velocidade

def calcula_parametros_caixa_particula(amplitude, num_onda):
    largura = 2 / math.pow(amplitude, 2)

    num_quantico = (num_onda * largura) / pi

    return largura, num_quantico

def calcula_comprimento_onda_broglie(velocidade, massa_particula):
    compr = cte_plank_joule / (massa_particula * velocidade)

    return compr

def converte_medida(medida, tipo_conversao):
    # converte de nanometros para metros
    medida_convert = medida * 1E-9

    # converte de metros para nanometros
    if tipo_conversao == 'nm':
        medida_convert = medida / 1E-9

    return medida_convert

def calculo_probab_unico_ponto(largura_caixa, num_quantico, ponto): 
    
    prob = (2 / largura_caixa) * (np.sin((num_quantico * pi / largura_caixa) * (ponto)) ** 2)

    return prob

def calculo_probab_dois_pontos(largura_caixa, num_quantico, ponto_a, ponto_b):
    x = np.linspace(ponto_a, ponto_b, 200)

    func_probab = (2 / largura_caixa) * (np.sin((num_quantico * pi / largura_caixa) * x) ** 2)

    probab = np.trapz(func_probab, x)

    return probab

def exibe_simulacao(largura_caixa, n_inicial, n_final, tipo_particula):
    energia_ni = calcula_energia_particula(largura_caixa, 1, tipo_particula)

    fig, ax = plt.subplots(figsize=(4,6))
    ax.set_xlim(0, largura_caixa)
    ax.set_xlabel('largura da caixa (m)')
    ax.set_ylim(0, math.pow(6, 2) * energia_ni)  # Definindo limite no eixo y
    ax.set_ylabel('Energia (J)')
    for nivel in range(1, 6):
        ax.plot([0, largura_caixa], [math.pow(nivel, 2) * energia_ni, math.pow(nivel, 2) * energia_ni], '-', color='C0', linewidth=1)
        fig.set_label(nivel)

    particula, = ax.plot([], [], 'bo', markersize=10)
    x = 0
    velocidade_x = largura_caixa / 50
    
    def init():
        particula.set_data([], [])
        return particula,

    def animate(frame):
        nonlocal x, velocidade_x
        x += velocidade_x
        
        
        if frame <= 100:
            y = math.pow(n_inicial, 2) * energia_ni
        else:
            y = math.pow(n_final, 2) * energia_ni            

        
        particula.set_data([x], [y])
        # Se a partícula atingir as extremidades do poço de potencial, muda a direção
        if x >= largura_caixa or x <= 0:
            velocidade_x *= -1
        return particula,


    ani = animation.FuncAnimation(fig, animate, frames=200, interval=25, blit=True, init_func=init)

    plt.show()

def main():
    imprime_tela_inicial()

    while True:
        menu_nav()

        opcao = int(input("\nopcao escolhida: "))

        if opcao == 1:
            tipo_particula = int(input("\nparticula a ser confinada:\n\t1 - eletron\n\t2 - proton\nopcao: "))

            largura_caixa = float(input("\ndigite a largura da caixa (L), em nanometros: "))
            
            n_inicial = int(input("digite o numero quantico inicial da particula (ni): "))
            n_final = int(input("digite o numero quantico final da particula (nf): "))
        
            largura_caixa = converte_medida(largura_caixa, 'm')

            amplitude_ni, num_onda_ni = calcula_parametros_funcao_onda(largura_caixa, n_inicial)
            amplitude_nf, num_onda_nf = calcula_parametros_funcao_onda(largura_caixa, n_final)

            print(f"\nfuncoes de onda:\n\tni: Ψ₁(x) = {amplitude_ni:.3G}sin({num_onda_ni:.3G}x)\n\tnf: Ψ₂(x) = {amplitude_nf:.3G}sin({num_onda_nf:.3G}x)")

            if tipo_particula == 1:
                energia_ni_joule = calcula_energia_particula(largura_caixa, n_inicial, 1)
                energia_nf_joule = calcula_energia_particula(largura_caixa, n_final, 1)

                print()

                print(f"energia da particula:\n\tni: {energia_ni_joule:.3G} J ou {(energia_ni_joule / 1.602e-19):.3G} eV\n\tnf: {energia_nf_joule:.3G} J ou {(energia_nf_joule / 1.602e-19):.3G} eV")

                energia_foton = calcula_energia_foton(n_inicial, n_final, 1, largura_caixa)
                tipo = "absorvido"

                if n_inicial > n_final:
                    tipo = "emitido"

                comp_onda, freq_foton = calcula_comprimento_onda_freq_foton(energia_foton)

                comp_onda_nm = converte_medida(comp_onda, 'nm')

                print()

                print(f"dados do foton {tipo.upper()}: \n\tenergia: {energia_foton:.3G} J ou {(energia_foton / 1.602e-19):.3G} eV\n\tcomprimento de onda: {comp_onda:.3G} m ou {comp_onda_nm:.3G} nm\n\tfrequencia: {freq_foton:.3G} Hz")

                vel_i = calcula_velocidade_particula(energia_ni_joule, massa_eletron)
                vel_f = calcula_velocidade_particula(energia_nf_joule, massa_eletron)

                print()

                print(f"velocidade da particula:\n\tni: {vel_i:.3G} m/s\n\tnf: {vel_f:.3G} m/s")

                compr_broglie_ni = calcula_comprimento_onda_broglie(vel_i, massa_eletron)
                compr_broglie_nf = calcula_comprimento_onda_broglie(vel_f, massa_eletron)

                compr_broglie_ni_conv = converte_medida(compr_broglie_ni, 'nm')
                compr_broglie_nf_conv = converte_medida(compr_broglie_nf, 'nm')

                print()

                print(f"comprimento de onda de De Broglie:\n\tni: {compr_broglie_ni:.3G} m ou {compr_broglie_ni_conv:.3G} nm\n\tnf: {compr_broglie_nf:.3G} m ou {compr_broglie_nf_conv:.3G} nm")

            elif tipo_particula == 2:
                energia_ni_joule = calcula_energia_particula(largura_caixa, n_inicial, 2)
                energia_nf_joule = calcula_energia_particula(largura_caixa, n_final, 2)

                print()

                print(f"energia da particula:\n\tni: {energia_ni_joule:.3G} J ou {(energia_ni_joule / 1.602e-19):.3G} eV\n\tnf: {energia_nf_joule:.3G} J ou {(energia_nf_joule / 1.602e-19):.3G} eV")

                energia_foton = calcula_energia_foton(n_inicial, n_final, 2, largura_caixa)
                tipo = "absorvido"

                if n_inicial > n_final:
                    tipo = "emitido"

                comp_onda, freq_foton = calcula_comprimento_onda_freq_foton(energia_foton)
                comp_onda_nm = converte_medida(comp_onda, 'nm')

                print()

                print(f"dados do foton {tipo.upper()}: \n\tenergia: {energia_foton:.3G} J.s ou {(energia_foton / 1.602e-19):.3G} eV\n\tcomprimento de onda: {comp_onda:.3G} m ou {comp_onda_nm:.3G} nm\n\tfrequencia: {freq_foton:.3G} Hz")

                vel_i = calcula_velocidade_particula(energia_ni_joule, massa_proton)
                vel_f = calcula_velocidade_particula(energia_nf_joule, massa_proton)

                print()

                print(f"velocidade da particula:\n\tni: {vel_i:.3G} m/s\n\tnf: {vel_f:.3G} m/s")

                compr_broglie_ni = calcula_comprimento_onda_broglie(vel_i, massa_proton)
                compr_broglie_nf = calcula_comprimento_onda_broglie(vel_f, massa_proton)

                compr_broglie_ni_conv = converte_medida(compr_broglie_ni, 'nm')
                compr_broglie_nf_conv = converte_medida(compr_broglie_nf, 'nm')

                print()

                print(f"comprimento de onda de De Broglie:\n\tni: {compr_broglie_ni:.3G} m ou {compr_broglie_ni_conv:.3G} nm\n\tnf: {compr_broglie_nf:.3G} m ou {compr_broglie_nf_conv:.3G} nm")

                
            exibe_simulacao(largura_caixa, n_inicial, n_final, tipo_particula)

            print("\ngostaria de calcular a probabilidade de encontrar a particula entre 2 pontos:\n\t1 - sim\n\t2 - nao")
            opcao = int(input("opcao: "))
            
            if opcao == 1:
                print("\ndigite as coordenadas de onde a particula sera procurada: ")
                
                while True:
                    coord_a = float(input("coordenada a, em nanometros: "))
                    coord_a = converte_medida(coord_a, 'm')
                    if coord_a <= 0 or coord_a > largura_caixa:
                        print("coordenada fora dos limites da caixa, por favor digite novamente.")
                    else:
                        break
                while True:
                    coord_b = float(input("coordenada b, em nanometros: "))
                    coord_b = converte_medida(coord_b, 'm')
                    if coord_b <= 0 or coord_b > largura_caixa:
                        print("coordenada fora dos limites da caixa, por favor digite novamente.")
                    else:
                        break

                probab_ni = calculo_probab_dois_pontos(largura_caixa, n_inicial, coord_a, coord_b)
                probab_nf = calculo_probab_dois_pontos(largura_caixa, n_final, coord_a, coord_b)

                print()

                print(f"probabilidade de encontrar a particula entre a e b:\n\tni: {(probab_ni * 100):.3G}%\n\tnf: {(probab_nf * 100):.3G}%")

            plota_graficos_func_onda(amplitude_ni, amplitude_nf, num_onda_ni, num_onda_nf, largura_caixa, n_inicial, n_final)

            plota_graficos_distr_probabilidade(amplitude_ni, amplitude_nf, num_onda_ni, num_onda_nf, largura_caixa, n_inicial, n_final)
        
        elif opcao == 2:
            print("\npara determinar parametros da caixa e da particula")

            input_amplitude = float(input("digite o valor da amplitude (A), em metros: "))
            input_num_onda = float(input("digite o valor do numero de onda (k): "))

            largura, num_quantico = calcula_parametros_caixa_particula(input_amplitude, input_num_onda)

            largura_conv_nm = converte_medida(largura, 'nm')

            print(f"\nparametros da caixa e da particula:\n\tlargura da caixa: {largura:.3G} m ou {largura_conv_nm:.3G} nm\n\tnivel quantico: {round(num_quantico)}")
            
            while True:
                posicao = float(input("\ndigite uma posicao [px]*L para determinar a probabilidade \nde a particula se encontrar nesta posicao: "))
                
                if (posicao * largura) > largura:
                    print("\nposicao localizada fora das dimensoes da caixa. favor digitar novamente")
                else:
                    break
            
            ponto_px = posicao * largura
            probab = calculo_probab_unico_ponto(largura, num_quantico, ponto_px)

            print(f"\nprobabilidade da particula estar em {posicao}*L: {probab:.3G}dx")

        elif opcao == 3:
            break
        else:
            print("\nopcao invalida! tente novamente.")
main()