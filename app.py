from flask import Flask, request, render_template_string, redirect, url_for
import math

app = Flask(__name__)

@app.errorhandler(404)
def page_not_found(e):
    return redirect(url_for('home'))

@app.route('/')
def home():
    return render_template_string('''
        <!DOCTYPE html>
        <html>
        <head>
            <title>Accueil</title>
            <style>
                body { 
                    font-family: Arial, sans-serif; 
                    text-align: center; 
                    padding: 50px; 
                    background-color: #f5f5f5;
                }
                .container {
                    max-width: 800px;
                    margin: 0 auto;
                    padding: 20px;
                }
                .welcome-box {
                    background-color: #ffffff;
                    padding: 40px;
                    border-radius: 15px;
                    box-shadow: 0 0 20px rgba(0,0,0,0.1);
                }
                .button {
                    background-color: #4CAF50;
                    color: white;
                    padding: 15px 30px;
                    text-decoration: none;
                    border-radius: 8px;
                    font-size: 18px;
                    transition: background-color 0.3s;
                }
                .button:hover {
                    background-color: #45a049;
                }
                .error-message {
                    color: #dc3545;
                    margin: 10px 0;
                    font-weight: bold;
                }
                .input-error {
                    border: 2px solid #dc3545 !important;
                }
                .form-group {
                    margin: 15px 0;
                }
                input {
                    padding: 8px;
                    border: 1px solid #ddd;
                    border-radius: 4px;
                    width: 200px;
                }
                .result-box {
                    background: #e9f7ef;
                    padding: 20px;
                    margin: 30px auto;
                    border-radius: 10px;
                    width: 80%;
                }
            </style>
        </head>
        <body>
            <div class="container">
                <div class="welcome-box">
                    <h1>Bonjour les PICs! 🧪</h1>
                    <p>Application de calcul des coefficients de diffusion UNIQUAC</p>
                    <a href="{{ url_for('diffusion_app') }}" class="button">Commencer</a>
                </div>
            </div>
        </body>
        </html>
    ''')

@app.route('/diffusion', methods=['GET', 'POST'])
def diffusion_app():
    default_values = {
        'x_A': '0.25',
        'D_AB_0': '2.1e-5',
        'D_BA_0': '2.67e-5',
        'q_A': '1.432',
        'q_B': '1.4',
        'r_A': '1.4311',
        'r_B': '0.92',
        'a_AB': '-10.7575',
        'a_BA': '194.5302',
        'T': '313.13',
        'D_exp': '1.33e-5'
    }

    if request.method == 'POST':
        error_messages = []
        error_fields = []
        form_data = request.form.to_dict()

        # Validation des entrées
        for key in default_values:
            value = form_data.get(key, '').strip()
            
            if not value:
                form_data[key] = default_values[key]
                continue
                
            try:
                cleaned_value = value.replace(',', '.')
                float(cleaned_value)
                form_data[key] = value  # Garde le format original
            except ValueError:
                error_messages.append(f"'{value}' n'est pas valide pour {key}")
                error_fields.append(key)
                form_data[key] = value

        if error_messages:
            return render_template_string('''
                <div class="container">
                    <h1>Calculateur de Diffusion</h1>
                    <form method="post">
                        {% for key in default.keys() %}
                        <div class="form-group">
                            <label>{{key}} :</label><br>
                            <input type="text" 
                                   name="{{key}}" 
                                   value="{{ form_data[key] }}"
                                   {% if key in error_fields %}class="input-error"{% endif %}>
                        </div>
                        {% endfor %}
                        <button type="submit" class="button">Calculer</button>
                    </form>
                    <div class="error-message">
                        {% for msg in error_messages %}
                        <p>{{ msg }}</p>
                        {% endfor %}
                        <div style="color:red; margin-top:10px;">Veuillez entrer uniquement des nombres valides !{{ error }}</div>
                    </div>
                    <a href="/" class="button">Retour à l'accueil</a>
                </div>
            ''', 
            default=default_values,
            form_data=form_data,
            error_messages=error_messages,
            error_fields=error_fields)

        # Conversion et calculs
        try:
            data = {key: float(form_data[key].replace(',', '.')) for key in default_values}
            
            # Calculs scientifiques
            x_B = 1 - data['x_A']
            y_A = data['r_A'] ** (1/3)
            y_B = data['r_B'] ** (1/3)
            phi_A = (data['x_A'] * y_A) / (data['x_A'] * y_A + x_B * y_B)
            phi_B = x_B * y_B / (data['x_A'] * y_A + x_B * y_B)
            theta_A = (data['x_A'] * data['q_A']) / (data['x_A'] * data['q_A'] + x_B * data['q_B'])
            theta_B = 1 - theta_A
            tau_AB = math.exp(-data['a_AB'] / data['T'])
            tau_BA = math.exp(-data['a_BA'] / data['T'])
            
            theta_AA_val = (theta_A * 1.0) / (theta_A * 1.0 + theta_B * tau_BA)
            theta_BB_val = (theta_B * 1.0) / (theta_A * tau_AB + theta_B * 1.0)
            theta_AB_val = (theta_A * tau_AB) / (theta_A * tau_AB + theta_B * 1.0)
            theta_BA_val = (theta_B * tau_BA) / (theta_A * 1.0 + theta_B * tau_BA)
            
            ln_D = (data['x_A'] * math.log(data['D_BA_0']) + x_B * math.log(data['D_AB_0']))
            ln_D += 2 * (data['x_A'] * math.log(data['x_A']/phi_A) + x_B * math.log(x_B/phi_B))
            ln_D += 2 * data['x_A'] * x_B * ((phi_A/data['x_A'])*(1 - y_A/y_B) + (phi_B/x_B)*(1 - y_B/y_A))
            ln_D += data['x_A'] * data['q_B'] * ((1 - theta_AB_val**2)*math.log(tau_AB) + (1 - theta_AA_val**2)*tau_BA*math.log(tau_BA))
            ln_D += x_B * data['q_A'] * ((1 - theta_BA_val**2)*math.log(tau_BA) + (1 - theta_BB_val**2)*tau_AB*math.log(tau_AB))
            
            D_calc = math.exp(ln_D)
            erreur = abs((D_calc - data['D_exp'])/data['D_exp']) * 100
            
            return render_template_string('''
                <div class="container">
                    <h1>Résultats du Calcul</h1>
                    <div class="result-box">
                        <h3>📊 Résultats :</h3>
                        <p>θ_AA = {{ "%.4f"|format(theta_AA) }}</p>
                        <p>θ_AB = {{ "%.4f"|format(theta_AB) }}</p>
                        <p>θ_BA = {{ "%.4f"|format(theta_BA) }}</p>
                        <p>θ_BB = {{ "%.4f"|format(theta_BB) }}</p>
                        <p>Coefficient calculé : {{ "%.2e"|format(D_calc) }} cm²/s</p>
                        <p>Valeur expérimentale : {{ "%.2e"|format(D_exp) }} cm²/s</p>
                        <p>Écart : {{ "%.2f"|format(erreur) }}%</p>
                    </div>
                    <a href="/diffusion" class="button">Nouveau Calcul</a>
                    <a href="/" class="button">Accueil</a>
                </div>
            ''',
            theta_AA=theta_AA_val,
            theta_AB=theta_AB_val,
            theta_BA=theta_BA_val,
            theta_BB=theta_BB_val,
            D_calc=D_calc,
            D_exp=data['D_exp'],
            erreur=erreur)

        except Exception as e:
            return render_template_string('''
                <div class="container">
                    <div class="error-message">
                        <h2>Erreur de Calcul ❌</h2>
                        <p>{{ error }}</p>
                    </div>
                    <a href="/diffusion" class="button">Réessayer</a>
                    <a href="/" class="button">Accueil</a>
                </div>
            ''', error=str(e))

    # GET Request
    return render_template_string('''
        <div class="container">
            <h1>Calculateur de Diffusion</h1>
            <form method="post">
                {% for key, value in default.items() %}
                <div class="form-group">
                    <label>{{key}} :</label><br>
                    <input type="text" name="{{key}}" value="{{ value }}">
                </div>
                {% endfor %}
                <button type="submit" class="button">Calculer</button>
            </form>
            <a href="/" class="button">Retour à l'accueil</a>
        </div>
    ''', default=default_values)

if __name__ == '__main__':
    app.run(debug=True)