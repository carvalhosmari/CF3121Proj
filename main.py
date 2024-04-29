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

def calcula_funcao_onda(largura_caixa, num_quantico):
    amplitude = math.sqrt(largura_caixa / 2)
    num_onda = (num_quantico * pi) / largura_caixa

    x = np.linspace(0, num_onda * largura_caixa, 100)

    funcao_onda = amplitude * np.sin(x)

    return funcao_onda, amplitude, num_onda

def calcula_energia_particula(largura_caixa, num_quantico, tipo_particula):
    energia_joule = math.pow(num_quantico, 2) * math.pow(cte_plank_joule, 2) / (8 * math.pow(largura_caixa, 2) * massa_eletron)

    # se particula confinada for proton
    if tipo_particula == 2:
        energia_joule = math.pow(num_quantico, 2) * math.pow(cte_plank_joule, 2) / (8 * math.pow(largura_caixa, 2) * massa_proton)        

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

def calcula_comprimento_onda_foton(energia_foton):
    comprimento_onda = (cte_plank_joule * velocidade_luz) / energia_foton

    return comprimento_onda

def calcula_velocidade_particula(energia_particula, massa_particula):
    velocidade = math.sqrt((2 * energia_particula) / massa_particula)

    return velocidade

def main():
    imprime_tela_inicial()
    
    while True: 
        menu_nav()

        opcao = int(input("\nopcao escolhida: "))

        if opcao == 1:
            n_inicial = int(input("digite o numero quantico inicial da particula (ni): "))
            n_final = int(input("digite o numero quantico final da particula (nf): "))
            largura_caixa = int(input("digite a largura da caixa (L), em nanometros: "))
            print("digite as coordenadas de onde a particula sera procurada: ")
            coord_a = int(input("coordenada a: "))
            coord_b = int(input("coordenada b: "))

        elif opcao == 2:
            # testes
            energia = calcula_energia_particula(2e-10, 1, 1)
            print(energia)
            compr = calcula_comprimento_onda_foton(energia)
            print(compr)
            velocidade = calcula_velocidade_particula(energia, 9.11e-31)
            print(velocidade)
        elif opcao == 3:
            break
        else:
            print("\nopcao invalida! tente novamente.")
main()