import math
import numpy as np
import matplotlib.pyplot as plt

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

    print("\n--------------------------------------------------------------------------------------------------------------")

def menu_nav():
    print("\nSelecione a entrada:\n\n\t1 - Determinar funcao de onda quantica e demais parametros;\n\t2 - Determinar parametros da particula e caixa a partir da funcao de onda\n\t3 - Encerrar programa")

def calcula_parametros_funcao_onda(largura_caixa, num_quantico):
    amplitude = math.sqrt((2) / (largura_caixa * 1E-9))
    num_onda = (num_quantico * pi) / (largura_caixa * 1E-9)

    return amplitude, num_onda

def plota_graficos_func_onda(amplitude_ni, amplitude_nf, num_onda_ni, num_onda_nf, largura_caixa, n_inicial, n_final):
    x = np.linspace(0, (largura_caixa * 1E-9), 150)

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
    
def calcula_energia_particula(largura_caixa, num_quantico, tipo_particula):
    energia_joule = (math.pow(num_quantico, 2) * math.pow(cte_plank_joule, 2)) / (8 * math.pow(largura_caixa * 1E-9, 2) * massa_eletron)

    # se particula confinada for proton
    if tipo_particula == 2:
        energia_joule = (math.pow(num_quantico, 2) * math.pow(cte_plank_joule, 2)) / (8 * math.pow(largura_caixa * 1E-9, 2) * massa_proton)        

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
def main():
    imprime_tela_inicial()
    
    while True: 
        menu_nav()

        opcao = int(input("\nopcao escolhida: "))

        if opcao == 1:
            n_inicial = int(input("digite o numero quantico inicial da particula (ni): "))
            n_final = int(input("digite o numero quantico final da particula (nf): "))
            largura_caixa = float(input("digite a largura da caixa (L), em nanometros: "))
            print("digite as coordenadas de onde a particula sera procurada: ")
            coord_a = int(input("coordenada a: "))
            coord_b = int(input("coordenada b: "))

            tipo_particula = int(input("particula a ser confinada:\n\t1 - eletron\n\t2 - proton\nopcao:"))

            amplitude_ni, num_onda_ni = calcula_parametros_funcao_onda(largura_caixa, n_inicial)
            amplitude_nf, num_onda_nf = calcula_parametros_funcao_onda(largura_caixa, n_final)

            plota_graficos_func_onda(amplitude_ni, amplitude_nf, num_onda_ni, num_onda_nf, largura_caixa, n_inicial, n_final)

            print(f"funcoes de onda:\n\tni: Ψ₁(x) = {amplitude_ni:.2E}sin({num_onda_ni:.2E}x)\n\tnf: Ψ₂(x) = {amplitude_nf:.2E}sin({num_onda_nf:.2E}x)")

            if tipo_particula == 1:
                energia_ni_joule = calcula_energia_particula(largura_caixa, n_inicial, 1)
                energia_nf_joule = calcula_energia_particula(largura_caixa, n_final, 1)
                
                print()

                print(f"energia da particula:\n\tni: {energia_ni_joule:.2E} J ou {(energia_ni_joule / 1.602e-19):.2f} eV\n\tnf: {energia_nf_joule:.2E} J ou {(energia_nf_joule / 1.602e-19):.2f} eV")
                
                energia_foton = calcula_energia_foton(n_inicial, n_final, 1, largura_caixa)
                tipo = "absorvido"

                if n_inicial > n_final:
                    tipo = "emitido"

                comp_onda, freq_foton = calcula_comprimento_onda_freq_foton(energia_foton)
                
                print()

                print(f"dados do foton {tipo.upper()}: \n\tenergia: {energia_foton:.2E} J ou {(energia_foton / 1.602e-19):.2f} eV\n\tcomprimento de onda: {comp_onda:.2E} m\n\tfrequencia: {freq_foton:.2E} Hz")

                vel_i = calcula_velocidade_particula(energia_ni_joule, massa_eletron)
                vel_f = calcula_velocidade_particula(energia_nf_joule, massa_eletron)

                print()

                print(f"velocidade da particula:\n\tni: {vel_i} m/s\n\tnf: {vel_f} m/s")

                compr_broglie_ni = calcula_comprimento_onda_broglie(vel_i, massa_eletron)
                compr_broglie_nf = calcula_comprimento_onda_broglie(vel_f, massa_eletron)

                print()

                print(f"comprimento de onda de De Broglie:\n\tni: {compr_broglie_ni:.2E} m\n\tnf:{compr_broglie_nf:.2E} m")
            elif tipo_particula == 2:
                energia_ni_joule = calcula_energia_particula(largura_caixa, n_inicial, 2)
                energia_nf_joule = calcula_energia_particula(largura_caixa, n_final, 2)
                
                print()

                print(f"energia da particula:\n\tni: {energia_ni_joule:.2E} J ou {(energia_ni_joule / 1.602e-19):.2f} eV\n\tnf: {energia_nf_joule:.2E} J ou {(energia_nf_joule / 1.602e-19):.2f} eV")
                
                energia_foton = calcula_energia_foton(n_inicial, n_final, 2, largura_caixa)
                tipo = "absorvido"

                if n_inicial > n_final:
                    tipo = "emitido"

                comp_onda, freq_foton = calcula_comprimento_onda_freq_foton(energia_foton)
                
                print()

                print(f"dados do foton {tipo.upper()}: \n\tenergia: {energia_foton:.2E} J.s ou {(energia_foton / 1.602e-19):.2f} eV\n\tcomprimento de onda: {comp_onda:.2E} m\n\tfrequencia: {freq_foton:.2E} Hz")
                
                vel_i = calcula_velocidade_particula(energia_ni_joule, massa_proton)
                vel_f = calcula_velocidade_particula(energia_nf_joule, massa_proton)

                print()

                print(f"velocidade da particula:\n\tni: {vel_i:.2E} m/s\n\tnf: {vel_f:.2E} m/s")

                compr_broglie_ni = calcula_comprimento_onda_broglie(vel_i, massa_proton)
                compr_broglie_nf = calcula_comprimento_onda_broglie(vel_f, massa_proton)

                print()

                print(f"comprimento de onda de De Broglie:\n\tni: {compr_broglie_ni:.2E} m\n\tnf: {compr_broglie_nf:.2E} m")

        elif opcao == 2:
            # testes
            largura, num_quantico = calcula_parametros_caixa_particula(5.612e4, 29.68e19)

            print(f"{largura:.2E} / {num_quantico}")
        elif opcao == 3:
            break
        else:
            print("\nopcao invalida! tente novamente.")
main()