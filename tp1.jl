### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ d1493290-93b6-451b-8974-eba1b68b7c0b
begin
	using LinearAlgebra, PlutoUI, Plots, DelimitedFiles
	TableOfContents()
end

# ╔═╡ 2532cd20-dd4c-11ee-2500-c37c88f25389
md"""
# TP 1: Ajuste de círculos
"""

# ╔═╡ c415e6b2-31b7-45c6-b6dc-eec5bdccd7d2
md"""
El objetivo de este trabajo es implementar y comparar algoritmos que ajusten conjuntos de datos con círculos, siguiendo diferentes metodologías. 
 A lo largo de todo el texto, asumiremos que los datos vienen dados por vectores ${\bf x} = (x_1,\dots,x_n)^t$ e ${\bf y} = (y_1,\dots,y_n)^t$, de manera tal que los puntos $(x_i,y_i)$ caigan aproximadamente en un círculo. 
"""

# ╔═╡ e6e60a47-3396-46ab-af1f-e4ab35216d46
md"""
## Ajuste algebraico

Un círculo puede definirse implíctamente con la ecuación:
```math
\alpha (x^2 + y^2) + \beta x + \gamma y + \delta = 0.
```
Para garantizar unicidad en la escritura, asumiremos $\alpha = 1$. Además, según veremos más adelante, resultará útil reescribir los parámetros, notando:
```math
x^2+y^2 -2ax-2by -c = 0.
```

Si todos los datos cayeran sobre un círculo tendríamos que la ecuación se satisface para todo $(x_i,y_i)$. Asumiendo que esto no sucede (los datos contienen error), buscaremos valores de $a$, $b$ y $c$ que minimicen:
```math
\sum_i (x_i^2 + y_i^2 -2ax_i-2by_i-c)^2.
```
O, equivalentemente:
```math
\min_{a,b,c} \big\|2a{\bf x} + 2b{\bf y} + c - {\bf p}\big\|, \quad\quad\textrm{donde }\; {\bf p} =\big(x_1^2+y_1^2, \dots,x_n^2+y_n^2\big)^t.
```

Este enfoque se conoce como *ajuste algebraico* del círculo, dado que se propone ajustar los parámetros de la ecuación algebraica que lo define.
"""

# ╔═╡ dc616049-f6d3-4478-983e-529dd5be1982
md"""
### Ejercicio 1

Los datos a los que le debemos calcular el ajuste (o con los cuales vamos a testear el ajuste) se encuentran en formato de archivo `.csv` dentro delarchivo `datos.zip`. Con la libreria `DelimitedFiles` podemos exportarlos para tenerlos en formato de matriz. A continuación se realiza la importación de uno de los conjuntos de datos como ejemplo. 

Importar los 10 conjuntos de datos y realizar un scatter plot de los mismos.
"""

# ╔═╡ 377661f1-0129-4e23-b82a-34fc0773b21c
begin
	# Definir la función para importar los datos y realizar el scatter plot
	function import_and_plot_data(file_name)
    	data = readdlm(file_name, ',')
    	x = data[:, 1]
	    y = data[:, 2]
    	scatter(x, y, label="", xlabel="X", ylabel="Y", title="Conjunto de Datos")
    	
    	# Extraer el número del archivo CSV
    	file_number = parse(Int, replace(split(split(file_name, "\\")[end], ".")[1], "datos" => ""))
		
	    # Guardar el gráfico como un archivo de imagen con el nombre adecuado
    	savefig("scatter_plots\\scatter_plot_$file_number.png")
	end
	
	# Lista de nombres de archivos de los conjuntos de datos con la nueva ruta
	file_names = ["C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos1.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos2.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos3.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos4.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos5.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos6.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos7.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos8.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos9.csv",
	              "C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos10.csv"]
	
	# Importar y plotear cada conjunto de datos
	for file_name in file_names
	    import_and_plot_data(file_name)
	end
end

# ╔═╡ 86f2a5d8-9901-47a3-81e2-9c5f351a25ae
md"""
### Ejercicio 2

Construir una función `def_problem` que tome como argumento a los datos y devuelva la matriz del ajuste junto con el vector del lado derecho teniendo en cuenta que lo que queremos es que estos se ajusten a una circunferencia.
"""

# ╔═╡ bc674ff7-e2ef-44f9-8fbe-673b05173177
begin
	function def_problem(data)
	    n = size(data, 1)
	    A = zeros(n, 3)
	    b = zeros(n)
	    
	    for i in 1:n
	        x_i = data[i, 1]
	        y_i = data[i, 2]
	        
	        A[i, 1] = 2 * x_i
	        A[i, 2] = 2 * y_i
	        A[i, 3] = 1
	        b[i] = x_i^2 + y_i^2
	    end
	    
	    return A, b
	end
	data = readdlm("C:\\Users\\Lenovo\\Desktop\\random\\uba\\optimizacion\\datos\\datos1.csv", ',')
	A, b12 = def_problem(data)

end

# ╔═╡ d5cc8b5c-b252-4587-a0f9-e3474159afbb
md"""
### Ejercicio 3

Implementar una función `ajuste_alg` que reciba como argumento los datos y resuelva el problema utilizando la descomposición `QR`. La salida deben ser los coeficientes $a$, $b$ y $c$.

"""

# ╔═╡ 8b01209a-2d94-4019-9983-dd40e745ab43
begin
	function ajuste_alg(data)
	    A, b = def_problem(data)
	    
	    # Aplicar descomposición QR reducida
	    Q0, R = qr(A)
		Q = Matrix{Float64}(Q0)
	    y = Q' * b
	    x = R \ y
	    
	    # Obtener los coeficientes a, b y c
	    a = x[1]
	    b = x[2]
	    c = x[3]
	    
	    return a, b, c
	end

	a, b, c = ajuste_alg(data)
	
end

# ╔═╡ 7c589215-5bfa-4cc3-a05d-c7bb6cb0c7e7
md"""
### Ejercicio 4

Implementar una función `plot_algebraico` que tome los coeficientes encontrados para un conjunto de datos y grafique el ajuste junto a los puntos del conjunto.
"""

# ╔═╡ 1e216b7d-0945-427e-a4fe-286e4e605b93
begin
	function plot_algebraico(a, b, c, data)
	    x = data[:, 1]
	    y = data[:, 2]
	    
	    # Calcular el centro y el radio del círculo ajustado
	    centro_x = a
	    centro_y = b
	    radio = sqrt(a^2 + b^2 + c)
	    
	    # Crear el círculo ajustado
	    θ = LinRange(0, 2π, 100)
	    circle_x = centro_x .+ radio * cos.(θ)
	    circle_y = centro_y .+ radio * sin.(θ)
	    
	    # Graficar el círculo ajustado y los puntos del conjunto de datos
	    scatter(x, y, label="Datos")
	    plot!(circle_x, circle_y, label="GN - Ajuste Algebraico", color="red")
	    xlabel!("X")
	    ylabel!("Y")
	    title!("Ajuste Algebraico de Círculo")
    	
		# Set aspect ratio to equal
    	plot!(aspect_ratio=:equal)
	end
	
	plot_algebraico(a, b, c, data)
end

# ╔═╡ 0770b04d-c4dd-4aaf-886e-ef451e876f92
md"""
### Ejercicio 5

Hacer los gráficos correspondientes a cada conjuntos de datos.
"""

# ╔═╡ 447093a3-41a5-415c-9a2d-404a369d7f74
begin

	# Realizar el ajuste y graficar para cada conjunto de datos
	for (i, file_name) in enumerate(file_names)
		# Importar los datos desde el archivo .csv
	    data = readdlm(file_name, ',')
	    
	    # Realizar el ajuste algebraico
	    a, b, c = ajuste_alg(data)

		# Definir la ruta de guardado del archivo
	    ruta_guardado = "plots_algebraicos/plot_algebraico$i.png"
	    
	    # Graficar el ajuste algebraico junto a los puntos del conjunto de datos
	    plot_algebraico(a, b, c, data)


		# Guardar el gráfico en la ruta especificada
	    savefig(ruta_guardado)
	end
	
end

# ╔═╡ 3640e8fd-1b6e-4c4a-9082-38b63e2bbf44
md"""
---
"""

# ╔═╡ c8df4623-2e32-4615-836f-9c94b801fcf9
md"""
## Ajuste geométrico
"""

# ╔═╡ 0f451e7e-6d37-4776-b8bf-23d3d4d98362
md"""
El ajuste algebraico, al ceñirse a un modelo de cuadrados mínimos clásico, puede resolverse de manera simple y veloz. Sin embargo, no tiene en cuenta la naturaleza geométrica del problema: no minimiza ninguna distancia. Esto hace que en algunos casos se obtengan resultados muy distintos de los esperados. El *ajuste geométrico* consiste en minimizar las distancias de los datos al círculo. 

 Notamos $C(a,b,r)$ al círculo con centro $(a,b)$ y radio $r$ y 
 ```math
d_i = \sqrt{(x_i-a)^2+(y_i-b)^2}-r,
```
 la distancia de $(x_i,y_i)$ a $C(a,b,r)$. El problema de determinar los parámetros $a,b,r$ que mejor ajusten los datos ${\bf x}$,${\bf y}$ está dado por:
```math
 \min F(a,b,r)
```
donde 

```math
F(a,b,r) = \sum_{i=1}^n d_i^2 = \sum_{i=1}^n\Big(\sqrt{(x_i-a)^2+(y_i-b)^2}-r\Big)^2.
```
Es importante observar que este es un problema de cuadrados mínimos *no lineal*, por lo cual no es posible aplicar de manera directa las técnicas habituales de cuadrados mínimos. 

Para resolver esto utilizaremos el método iterativo de Gauss-Newton, que describimos a continuación.
El objetivo es minimizar 
```math
\sum_{i=1}^n d_i(\theta)^2, \quad \theta = (a,b,r).
```
Notamos ${\bf d}(\theta) = (d_1(\theta),\dots,d_n(\theta))^t$.
Idealmente (si los datos coincidieran en un mismo círculo) tendríamos que la solución optima $\theta^*$ satisface ${\bf d}(\theta^*) = {\bf 0}$.

Desarrollando el polinomio de Taylor a orden 1 tenemos que:
```math
{\bf d}(\theta+{\bf h}) \simeq {\bf d}(\theta)+J(\theta){\bf h}\sim {\bf 0},
```
donde $J$ es la matriz diferencial de ${\bf d}$. El método de Gauss-Newton toma una solución inicial $\theta_0$, resuelve el problema de cuadrados mínimos:
```math
J(\theta_0){\bf h} = -{\bf d}(\theta_0), 
```
toma $\theta_1 = \theta_0 + {\bf h}$, e itera el procedimiento, calculando una sucesión $\theta_i$ y deteniéndose cuando se satisface algún criterio de convergencia. 
"""

# ╔═╡ 217120b1-8e05-4348-9a99-fde8f5be7efd
md"""
### Ejercicio 6

Implementar las siguientes funciones:
1. una función que retorne el vector ${\bf d}(\theta)$, 
2. una función que retorne el valor de $F$
3. una función que calcule la matriz diferencial $J$ de ${\bf d}$ dependiendo del $\theta_0$.
"""

# ╔═╡ d0e96c17-368d-4daf-a9f8-4fa58e85778c
begin
	function distancia(data, theta)
	    a, b, r = theta
	    n = size(data, 1)
	    d = zeros(n)
	    
	    for i in 1:n
	        x_i = data[i, 1]
	        y_i = data[i, 2]
	        d[i] = sqrt((x_i - a)^2 + (y_i - b)^2) - r
	    end
	    
	    return d
	end
	function valor_F(data, theta)
	    d = distancia(data, theta)
	    return sum(d .^ 2)
	end
	function matriz_J(data, theta0)
	    a0, b0, r0 = theta0
	    n = size(data, 1)
	    J = zeros(n, 3)
	    
	    for i in 1:n
	        x_i = data[i, 1]
	        y_i = data[i, 2]
	        
	        # Derivadas parciales de d respecto a a, b y r
	        d_da = -(x_i - a0) / sqrt((x_i - a0)^2 + (y_i - b0)^2)
	        d_db = -(y_i - b0) / sqrt((x_i - a0)^2 + (y_i - b0)^2)
	        d_dr = -1
	        
	        # Llenar la fila i de la matriz J
	        J[i, 1] = d_da
	        J[i, 2] = d_db
	        J[i, 3] = d_dr
	    end
	    
	    return J
	end
	
end

# ╔═╡ 9e2ef2ba-3414-4a42-b340-83462378f19c
md"""
### Ejercicio 7

Implementar una función `ajuste_geom` que implemente el algoritmo de Gauss-Newton descripto, resolviendo 
```math
J(\theta_0){\bf h} = -{\bf d}(\theta_0), 
```
a través de la descomposición en valores singulares de $J$ ($J=U\Sigma V^t$). La función deberá impirmir también un mensaje si el problema está mal condicionado, es decir si $J$ tiene un número de condición muy grande.

**Sugerencia:** El proceso iterativo debe comenzar con algun $\theta_0$. Este puede tomarlo como la solución del algebraico.
"""

# ╔═╡ b22f936a-9f46-44e9-814a-70fbf334ef8d
function ajuste_geom(data, theta0, max_iter = 100, tol = 1e-6, λ = 1e-3)
    iter = 0
    theta = theta0
    while true
        iter += 1
        
        # Calcular la matriz diferencial J y el vector de distancias d
        J = matriz_J(data, theta)
        d = distancia(data, theta)

        
        # Calcular la descomposición en valores singulares de J
		decomposicion = svd(J)
		U = decomposicion.U
		Σ = decomposicion.S
		Vt = decomposicion.Vt
        
        # Compute h
        h = Vt' * ((1 ./ Σ) .* (U' * (-d)))
        
        # Check for NaN values in h
        if any(isnan, h)
            println("NaN values encountered. Stopping iteration.")
            return theta
        end

        # Actualizar los parámetros a, b, r
        theta = (theta[1] + h[1], theta[2] + h[2], theta[3] + h[3])
        
        # Verificar la convergencia
        if norm(h) < tol || iter >= max_iter
            break
        end
    end
    return theta
end


# ╔═╡ 889aa225-8bd9-4706-9c49-6418f4bc4ef1
md"""
### Ejercicio 8

Realizar un plot del ajuste geometrico obtenido junto con el ajuste algebraico y los datos. ¿Observa alguna mejora?
"""

# ╔═╡ 81719f65-2c80-4b07-8cc5-079770d4b034
begin
  
	# Realizar el ajuste y graficar para cada conjunto de datos
	for (i, file_name) in enumerate(file_names)
		# Importar los datos desde el archivo .csv
	    data = readdlm(file_name, ',')
	    
	    # Realizar el ajuste algebraico
	    a, b, c = ajuste_alg(data)

		# Definir la ruta de guardado del archivo
	    ruta_guardado = "plots_algebraicos_ygeom/plot_algebraicos_ygeom$i.png"
	    
	    # Graficar el ajuste algebraico junto a los puntos del conjunto de datos
	    plot_algebraico(a, b, c, data)

		# Graficar el ajuste geometrico junto a los puntos del conjunto de datos
		a1, b1, radio = ajuste_geom(data, (a,b,c+a*a+b*b))

    	
		# Calcular el centro y el radio del círculo ajustado
	    centro_x = a1
	    centro_y = b1
	    
		# Crear el círculo ajustado
	    θ = LinRange(0, 2π, 100)
	    circle_x = centro_x .+ radio * cos.(θ)
	    circle_y = centro_y .+ radio * sin.(θ)
	
	    plot!(circle_x, circle_y, label="Ajuste Geometrico", color="blue")


		# Guardar el gráfico en la ruta especificada
	    savefig(ruta_guardado)
	end
	
end

# ╔═╡ 15b76fbf-5861-4f92-96c1-c3a84a61c810
# El ajuste algebraico asume que los datos son solo una parte de la poblacion y que la poblacion entera esta uniformemente distribuida en el contorno del circulo.
# EL ajuste geometrico asume que todos los datos que se dan son representativos y por lo tanto ajusta el circulo para que este lo mas cerca posible de todos los puntos posible. 

# ╔═╡ 9bd5c0e1-9a3a-461d-ad77-fd633ce728eb
md"""
## Datos singulares
El método de Gauss-Newton asume que $J$ es definida positiva, lo que implica, en particular, que es inversible. Sin embargo, esto no tiene por qué ser cierto (depende de los datos). Más aún, podría ocurrir que algún valor singular sea muy cercano a cero, lo que daría una $\Sigma$ muy mal condicionada, ocasionando graves errores numéricos. Es importante remarcar que este inconveniente no proviene del uso de la descomposición en valores singulares, sino que es intrínseco a la matriz $J$ y, más esencialmente, a los datos ${\bf x}$ e ${\bf y}$.  

Para evitar esto, proponemos dos alternativas, de las cuales debe elegirse sólo una:
"""

# ╔═╡ 0f28b4da-a147-4e1d-91e7-18c425ceedec
md"""
#### Levenberg-Marquardt


Levenberg y Marquardt propusieron una modificación del algoritmo que consiste reemplazar $\Sigma$ por
```math
\tilde{\Sigma} = \Sigma + \lambda I.
```
eligiendo $\lambda>0$ lo más chico posible de manera que el mínimo valor singular de $\tilde{\Sigma}$ resulte razonablemente lejano a cero.  Es decir: $\lambda = \sigma_{\min}+\delta$, donde $\sigma_{min}$ es el mínimo valor singular y $\delta$ es una constante positiva. 

Si $\sigma_{min}$ es positivo pero muy pequeño, una posibilidad de correción es dejar el mínimo valor singular en $2\sigma_{min}+\delta$. En este caso, $\tilde{\Sigma}$ tendrá un número de condición aceptable. Si $\sigma_{min}$ es positivo y no-chico, no es necesario corregir $\Sigma$. 

Vale la pena observar que cuanto mayor es el valor de $\lambda$, más se parece el vector ${\bf h}$ resultante al dado por el método del gradiente; cuanto menor es $\lambda$, más parecida es la iteración a la correspondiente al método de Gauss-Newton puro. 

_Sugerencia_: Tomar un valor inicial $\delta$ (por ejemplo $\delta = 10^{-8}$). Se considera la solución inicial $\theta=(a,b,r)$ y se calcula un valor $\theta'=(a',b',r')$ , correspondiente a realizar una iteración del método con la corrección. Para estas soluciones se computa el valor del funcional $F$. Si $\theta'$ arroja valores de $F$ mayores que el correspondiente a $\theta$, entonces se agranda $\lambda$ con algún criterio hasta hallar un $\lambda$ que mejore el valor del funcional. 

"""

# ╔═╡ 90d3d702-84c2-45df-ae21-b04cad4c433c
md"""
#### Factorización Positiva

A partir de un cuestionamiento al método anterior, se derivó la siguiente propuesta alternativa. En lugar de corregir  _todos_ los valores singulares, se corrigen sólo los que traen inconvenientes, pero de manera distinta. La idea es construir directamente $\tilde{\Sigma}^{-1}$, del siguiente modo. Llamemos $\sigma_i$ al $i-$ésimo valor singular. 

+ Si $|\sigma_i|<\delta$ para algún $\delta$ chico, entonces $(\tilde{\Sigma}^{-1})_{ii} = \frac{\delta}{|\sigma_{max}|}$.
+ En caso contrario: $(\tilde{\Sigma}^{-1})_{ii} = \frac{1}{|\sigma_i|}$.

Un valor razonablemente pequeño de $\delta$ puede ser $\delta = \sqrt{\varepsilon} \approx 10^{-8}$.
"""

# ╔═╡ a3453611-9ced-47a1-b2c5-5163f558405a
md"""
### Ejercicio 9

Implementar una función `ajuste_mejorado` que ejecute el método de Gauss-Newton con la variante de Levenberg-Marquardt o con el método de Factorización Positiva, descriptos anteriormente.

La función debe devolver la solución, pero también el número de pasos realizados (es decir, cuántas iteraciones fueron necesarias para alcanzar la solución).

"""

# ╔═╡ fe65d66b-006b-477f-b147-9fa36c172d11
begin

	function paso_levmarq(J, f, lambda)
	    m, n = size(J)
	    A = J' * J
	    b = J' * f
	    delta_theta = -(A + lambda * I) \ b
	    return delta_theta
	end
	
	function paso_fact_positiva(J, f, delta)
	    m, n = size(J)
	    U, Σ, V = svd(J)
	    Σp = [(σ < delta) ? delta : σ for σ in Σ]
	    delta_theta = V * Diagonal(1 ./ Σp) * U' * f
	    return delta_theta
	end
	
	function ajuste_mejorado(data, theta0, metodo, max_iter=1000, tol = 1e-6, delta=1e-8)
		iter = 0
		theta = theta0
		
	    while true	      
			iter += 1

	        # Calcular el jacobiano
	        J = matriz_J(data, theta)
		  	d = distancia(data, theta)
			
			# Calcular la descomposición en valores singulares de J
			decomposicion = svd(J)
			U = decomposicion.U
			Σ = decomposicion.S
			Vt = decomposicion.Vt

	        
	        # Calcular el paso según el método especificado
	        if metodo == "Levenberg-Marquardt"
				sigma_min = last(svdvals(Σ))
				lambda = delta + sigma_min
				Σ_reemplazo = Diagonal(Σ) + lambda*I
				
				# Compute h
				h = Vt' * ((1 ./ Σ) .* (U' * (-d)))
	        elseif metodo == "Factorización Positiva"
				m, n = size(J)
				U, Σ, V = svd(J)
				sigma_max = first(svdvals(Σ))

				Σp_inv = [(σ < delta) ? delta/sigma_max : 1/σ for σ in Σ]
				# Compute h
				h = Vt' * (Diagonal(Σp_inv) * (U' * (-d)))

	        else
	            error("Método no válido.")
	        end
			
			theta = (theta[1] + h[1], theta[2] + h[2], theta[3] + h[3])

		  	# Verificar la convergencia
			if norm(h) < tol || iter >= max_iter
			  break
			end
		end
	    
	    # Devolver los parámetros ajustados y el número de iteraciones
	    return theta, iter
	end

	
end

# ╔═╡ 9e09614b-56db-4a69-8198-565eeeeb89dc
md"""
### Ejercicio 10

Comparar gráficamente los resultados obtenidos entre alguno de los últimos métodos descriptos y el metodo de Gauss-Newton puro para los casos donde la matriz $J$ resultante está mal condicionada.
"""

# ╔═╡ caf062a8-ec85-4512-b104-62ababead6c7
begin

	
	# Realizar el ajuste y graficar para cada conjunto de datos
	for (i, file_name) in enumerate(file_names)
		# Importar los datos desde el archivo .csv
	    data = readdlm(file_name, ',')
	    
	    # Realizar el ajuste algebraico
	    a, b, c = ajuste_alg(data)

		# Graficar el ajuste geometrico junto a los puntos del conjunto de datos
		a_gn, b_gn, radio_gn = ajuste_geom(data, (a,b,c))

        # Realizar el ajuste mejorado con Levenberg-Marquardt
        theta_lm, iter_lm = ajuste_mejorado(data, (a_gn, b_gn, radio_gn), "Levenberg-Marquardt")
		
		a_lm, b_lm, c_lm = theta_lm

        # Realizar el ajuste mejorado con Factorización Positiva
        theta_fp, iter_fp = ajuste_mejorado(data, (a_gn, b_gn, radio_gn), "Factorización Positiva")
		a_fp, b_fp, c_fp = theta_fp
		
		# Definir la ruta de guardado del archivo
	    ruta_guardado = "plots_comparacion/plot_comparacion$i.png"

    	
		# Crear el círculo ajustado
	    θ = LinRange(0, 2π, 100)

		x = data[:, 1]
	    y = data[:, 2]

	    # Graficar el círculo ajustado y los puntos del conjunto de datos
	    scatter!(x, y, label="Datos")
	    
	    # Graficar el ajuste algebraico junto a los puntos del conjunto de datos
	    plot_algebraico(a, b, c, data)

        # Graficar el ajuste mejorado con Levenberg-Marquardt
        plot!(a_lm .+ c_lm * cos.(θ), b_lm .+ c_lm * sin.(θ), label="Levenberg-Marquardt_$iter_lm", color="blue")

        # Graficar el ajuste mejorado con Factorización Positiva
        plot!(a_fp .+ c_fp * cos.(θ), b_fp .+ c_fp * sin.(θ), label="Factorización Positiva_$iter_fp", color="green")

		# Set aspect ratio to equal
    	plot!(aspect_ratio=:equal)

		# Guardar el gráfico en la ruta especificada
	    savefig(ruta_guardado)
	end
	
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DelimitedFiles = "~1.9.1"
Plots = "~1.40.2"
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "364884a1cfda4e8855215f24cc190c1acfd113bc"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a4c43f59baa34011e303e76f5c8c91bf58415aaf"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "87944e19ea747808b73178ce5ebb74081fdf2d35"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "0f4b5d62a88d8f59003e43c25a8a90de9eb76317"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.18"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Format]]
git-tree-sha1 = "f3cf88025f6d03c194d73f5d13fee9004a108329"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.6"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "3437ade7073682993e092ca570ad68a2aba26983"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a96d5c713e6aa28c242b0d25c1347e258d6541ab"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.3+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "359a1ba2e320790ddbe4ee8b4d54a305c0ea2aff"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.0+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "db864f2d91f68a5912937af80327d288ea1f3aee"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cad560042a7cc108f5a4c24ea1431a9221f22c1b"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.2"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dae976433497a2f841baadea93d27e68f1a12a97"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.39.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0a04a1318df1bf510beb2562cf90fb0c386f58c4"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "af81a32750ebc831ee28bdaaba6e1067decef51e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "3c403c6590dd93b36752634115e20137e79ab4df"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.2"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "3caa21522e7efac1ba21834a03734c57b4611c7e"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.4"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "07e470dabc5a6a4254ffebc29a1b3fc01464e105"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "31c421e5516a6248dfb22c194519e37effbf1f30"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─2532cd20-dd4c-11ee-2500-c37c88f25389
# ╟─c415e6b2-31b7-45c6-b6dc-eec5bdccd7d2
# ╟─e6e60a47-3396-46ab-af1f-e4ab35216d46
# ╠═d1493290-93b6-451b-8974-eba1b68b7c0b
# ╟─dc616049-f6d3-4478-983e-529dd5be1982
# ╠═377661f1-0129-4e23-b82a-34fc0773b21c
# ╟─86f2a5d8-9901-47a3-81e2-9c5f351a25ae
# ╠═bc674ff7-e2ef-44f9-8fbe-673b05173177
# ╟─d5cc8b5c-b252-4587-a0f9-e3474159afbb
# ╠═8b01209a-2d94-4019-9983-dd40e745ab43
# ╟─7c589215-5bfa-4cc3-a05d-c7bb6cb0c7e7
# ╠═1e216b7d-0945-427e-a4fe-286e4e605b93
# ╟─0770b04d-c4dd-4aaf-886e-ef451e876f92
# ╠═447093a3-41a5-415c-9a2d-404a369d7f74
# ╟─3640e8fd-1b6e-4c4a-9082-38b63e2bbf44
# ╟─c8df4623-2e32-4615-836f-9c94b801fcf9
# ╟─0f451e7e-6d37-4776-b8bf-23d3d4d98362
# ╟─217120b1-8e05-4348-9a99-fde8f5be7efd
# ╠═d0e96c17-368d-4daf-a9f8-4fa58e85778c
# ╟─9e2ef2ba-3414-4a42-b340-83462378f19c
# ╠═b22f936a-9f46-44e9-814a-70fbf334ef8d
# ╟─889aa225-8bd9-4706-9c49-6418f4bc4ef1
# ╠═81719f65-2c80-4b07-8cc5-079770d4b034
# ╠═15b76fbf-5861-4f92-96c1-c3a84a61c810
# ╟─9bd5c0e1-9a3a-461d-ad77-fd633ce728eb
# ╟─0f28b4da-a147-4e1d-91e7-18c425ceedec
# ╟─90d3d702-84c2-45df-ae21-b04cad4c433c
# ╟─a3453611-9ced-47a1-b2c5-5163f558405a
# ╠═fe65d66b-006b-477f-b147-9fa36c172d11
# ╟─9e09614b-56db-4a69-8198-565eeeeb89dc
# ╠═caf062a8-ec85-4512-b104-62ababead6c7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
