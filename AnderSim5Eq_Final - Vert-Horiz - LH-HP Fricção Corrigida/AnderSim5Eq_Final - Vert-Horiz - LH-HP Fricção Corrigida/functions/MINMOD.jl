function MINMOD(Q1, Q2, Q3)
    som_sign = sign(Q1) + sign(Q2) + sign(Q3)
    som_sign = abs(som_sign)      # dará 3 ou 1

    resultado = sign(Q1)*minimum([abs(Q1),abs(Q2),abs(Q3)])*(som_sign-1)/2
    # só dará diferente de zero se os três argumentos tiverem o mesmo sinal

     return (resultado)
end