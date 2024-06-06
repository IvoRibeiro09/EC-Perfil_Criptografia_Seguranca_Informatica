def eval_vertical1(B, A):
    x_B, y_B = B[0], B[1]
    x_A, y_A = A[0], A[1]
    
    # Cálculo da linha vertical
    r = (x_B - x_A) 
    
    return r

def eval_tangent1(B, A):
    x_B, y_B = B[0], B[1]
    x_A, y_A = A[0], A[1]
    
    # Casos especiais
    if A == (0, 0):
        return (1, 0)  # Se A for o ponto de identidade, retorna 1
    
    if y_A == 0:
        # Se y_A for zero, retorna o resultado de EvalVertical1
        return eval_vertical1(B, A)
    
    # Cálculo dos coeficientes da linha
    a = (-3 * x_A**2) 
    b = (2 * y_A) 
    c = (-b * y_A - a * x_A) 
    
    # Avaliação em B
    r = (a * x_B + b * y_B + c) 
    
    return r

def eval_line1(B, A_prime, A_double_prime):
    x_B, y_B = B[0], B[1]
    x_prime, y_prime = A_prime
    x_double_prime, y_double_prime = A_double_prime
    
    # Casos especiais
    if A_prime == (0, 0):
        return eval_vertical1(B, A_double_prime)
    
    if A_double_prime == (0, 0):
        return eval_vertical1(B, A_prime)
    
    if A_prime == (-x_double_prime, -y_double_prime):
        return eval_vertical1(B, A_prime)
    
    if A_prime == A_double_prime:
        return eval_tangent1(B, A_prime)
    
    # Cálculo dos coeficientes da linha
    a = (y_prime - y_double_prime) 
    b = (x_double_prime - x_prime) 
    c = (-b * y_prime - a * x_prime) 
    
    # Avaliação em B
    r = (a * x_B + b * y_B + c) 
    
    return r

def projective_point_double1(A):
    x, y, z = A[0], A[1], A[2]
    
    # Caso especial: ponto de identidade
    if z == 0 or y == 0:
        return (0, 1, 0)
    
    # Cálculo das coordenadas projetivas dobradas
    lambda_1 = (3 * x**2) 
    z_prime = (2 * y * z) 
    lambda_2 = (y**2) 
    lambda_3 = (4 * lambda_2 * x) 
    x_prime = (lambda_1**2 - 2 * lambda_3) 
    lambda_4 = (8 * lambda_2**2) 
    y_prime = (lambda_1 * (lambda_3 - x_prime) - lambda_4) 
    
    return (x_prime, y_prime, z_prime)



def pairing(E, p, q, A, B):
    if E.is_supersingular():
        return pairing1(E, p, q, A, B)
    else:
        raise ValueError("Only type-1 curves are supported")

def pairing1(E, p, q, A, B):
    x, y = B[0], B[1]
    
    a_zeta = (p - 1) // 2
    b_zeta = pow(3, (p + 1) // 4, p)
    
    x_prime = x * a_zeta % p
    
    B_prime = (x_prime, y)
    
    
    return tate_pairing(E,p,q,A,B_prime)

def tate_pairing(E,p,q,A,B):
    if E.is_supersingular():
        return tate_miller_solinas(A, B, p, q)
    else:
        raise ValueError("Only type-1 curves are supported")


def tate_miller_solinas(A, B, p, q):
    a, b, s, c = decompose_q(q)
    # Inicialização
    v_num = 1
    v_den = 1
    V = (A[0], A[1], 1)  
    t_num = 1
    t_den = 1
    
    # Contribuição (s * 2^b)
    for n in range(b):
        t_num = t_num**2
        t_den = t_den**2
    #verificar os proximos 3
        t_num = t_num * eval_tangent1(B, (V[0]/ V[2]**2, V[1]/ V[2]**3))
        V = projective_point_double1(V)
        t_den = t_den * eval_vertical1(B, (V[0]/ V[2]**2, V[1]/ V[2]**3))
    
    # Normalização
    V_b = (V[0] / V[2]**2) , (s * V[1] / V[2]**3)
    
    # Acumulação
    if s == -1:
        v_num = v_num * t_den
        v_den = v_den * t_num * eval_vertical1(B, (V[0]/ V[2]**2, V[1]/ V[2]**3))
    elif s == 1:
        v_num = v_num * t_num
        v_den = v_den * t_den
    
    # Contribuição 2^a
    for n in range(b, a -1):
        t_num = t_num**2
        t_den = t_den**2
        t_num = t_num * eval_tangent1(B, (V[0] / V[2]**2, V[1]/ V[2]**3))
        V = 2 *  projective_point_double1(V)
        t_den = t_den * eval_vertical1(B, (V[0] / V[2]**2, V[1]/ V[2]**3))
    
    # Normalização
    V_a = V[0] / V[2]**2, s * V[0] / V[2]**3
    
    # Acumulação
    v_num = v_num * t_num
    v_den = v_den * t_den
    
    # Correção para as contribuições (s * 2^b) e (c)
    v_num = v_num * eval_line1(B, V_a, V_b)
    v_den = v_den * eval_vertical1(B,(V_a[0]+V_b[0],V_a[1]+V_b[1]))
    if c == -1:
        v_den = v_den* eval_vertical1(B, A)
    
    # Correção do expoente
    eta = (p**2 - 1)/ q
    
    # Resultado final
    return (v_num / v_den)**eta

def decompose_q(q):
    a = 0
    b = 0
    s = 0
    c = 0
    
    # Encontra o valor de a
    while q % 2 == 0:
        a += 1
        q //= 2
    
    # Encontra o valor de s
    if q != 1:
        s = 1
        b = 0
        while (q - 1) % 2 == 0:
            b += 1
            q = (q - 1) // 2
        c = (q - 1) // 2
        
        # Ajusta c se necessário
        if c > 1:
            c = c - 2
        elif c == 0:
            c = -1
    
    return a, b, s, c

import hashlib
from math import ceil

def Hashfunc(hashfcn):
    oid_map = {
        "1.3.14.3.2.26": hashlib.sha1,
        "2.16.840.1.101.3.4.2.4": hashlib.sha224,
        "2.16.840.1.101.3.4.2.1": hashlib.sha256,
        "2.16.840.1.101.3.4.2.2": hashlib.sha384,
        "2.16.840.1.101.3.4.2.3": hashlib.sha512
    }
    return oid_map[hashfcn]
##direita
def HashToPoint(E, p, q, id, hashfcn):
    if E.is_supersingular() or (not q.is_prime()) or (not p.is_prime()):
        return HashToPoint1(E,p, q, id, hashfcn)
    else:
        raise ValueError("Only type-1 curves are supported")
##direita
def HashToPoint1(E,p, q, id, hashfcn):
    # Step 1: Hash the identity string id to an element of F_p
    y = HashToRange(id.encode(), p, hashfcn)
    y_square = pow(y,2,p)-1
    exp = (2 * p - 1) // 3
    # Step 2: Calculate x = (y^2 - 1)^((2 * p - 1) / 3) modulo p
    x =  pow(y_square,exp,p)
    
    # Step 3: Create the point Q' = (x, y)
    if x==0 and y ==0:
        raise ValueError("Q' most be a non-zero point in E(F_p)")
    Q_prime = E(x, y)
    
    # Step 4: Calculate Q = [(p + 1) / q] * Q', a point of order q in E(F_p)
    scalar = (p + 1) / q

    Q = scalar * Q_prime

    return Q

###direita
def HashToRange(s, n, hashfcn):
    # Get the length of the output of hash function hashfcn
    hashlen = hashfcn().digest_size
    
    v_ = [None] * 3
    h_ = [None] * 3
    # Initialize v_0 and h_0
    v_[0] = 0
    h_[0] = b'\x00' * hashlen
    
    # Iterate for i = 1 to 2
    for i in range(1, 3):
        # Concatenate h_(i - 1) and s
        t_i = h_[i-1] + s
        
        # Calculate h_i = hashfcn(t_i)
        h_[i] = hashfcn(t_i).digest()
        
        # Convert h_i to an integer a_i in the range 0 to 256^hashlen - 1
        a_i = int.from_bytes(h_[i], byteorder='big')
        
        # Calculate v_i = 256^hashlen * v_(i - 1) + a_i
        v_[i] = 256**hashlen * v_[i-1] + a_i
    
    # Calculate v = v_l (mod n)
    v = v_[2] % n
    
    return v

def HashBytes(b, p, hashfcn):
    # Step 1: Determine the length of the hash output
    hashlen = hashfcn().digest_size

    # Step 2: Compute the hash of the key string p
    K = hashfcn(p).digest()

    # Step 3: Initialize h_0 with null octets
    h_0 = b'\x00' * hashlen

    # Step 4: Calculate the number of iterations needed
    l = ceil(b / hashlen)  # Ceiling division

    # Step 5: Iterate l times
    r = b''
    for i in range(1, l + 1):
        # Step 5(a): Compute h_i = hashfcn(h_(i-1))
        h_i = hashfcn(h_0 if i == 1 else h_i).digest()

        # Step 5(b): Compute r_i = hashfcn(h_i || K)
        r_i = hashfcn(h_i + K).digest()

        # Append r_i to r
        r += r_i

    # Step 6: Take the leftmost b octets of r
    result = r[:b]

    return result

from math import ceil, log2
log2(4)
def canonical_encoding(E,p, k, o,v):
    if E.is_supersingular():
        return canonical_encoding1(v, p, k, o)
    else:
        raise ValueError("Only type-1 curves are supported")
        
def canonical_encoding1(v, p, k, o):
    l = ceil(log2(p) / 8)  # Ceiling(lg(p) / 8)

    # Step 2: Split v into real and imaginary parts
    i = pow(-1, (p - 1) // 4, p)
    b = v // i
    a = v - b * i
    if int(a) + b * i != v:
        raise ValueError("a + b * i must be equals v")
    
    # Step 3: Encode the real part (a) as a big-endian zero-padded fixed-length octet string
    a_256l = int(a).to_bytes(l, byteorder='big', signed=False)

    # Step 4: Encode the imaginary part (b) as a big-endian zero-padded fixed-length octet string
    b_256l = int(b).to_bytes(l, byteorder='big', signed=False)

    # Step 5: Depending on the choice of ordering o
    if o == 0:
        # If o = 0, concatenate a_(256^l) followed by b_(256^l)
        s = a_256l + b_256l
    elif o == 1:
        # If o = 1, concatenate b_(256^l) followed by a_(256^l)
        s = b_256l + a_256l
    else:
        raise ValueError("Invalid ordering parameter. It should be either 0 or 1")

    # Step 6: Return the encoded string
    return s