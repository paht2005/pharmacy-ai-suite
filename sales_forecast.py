# === Sales Forecast (Prophet) ===
# Phat Nguyen Cong - Github: https://github.com/paht2005

import pandas as pd
from prophet import Prophet
import matplotlib.pyplot as plt

# Load or simulate sales data (for testing)
def load_sales_data():
    data = {
        'ds': pd.date_range(start="2024-01-01", periods=120),
        'y': [200 + i % 30 - (i // 10 % 3)*5 for i in range(120)]
    }
    return pd.DataFrame(data)

def forecast_sales(df, future_days=30):
    model = Prophet(daily_seasonality=True)
    model.fit(df)
    future = model.make_future_dataframe(periods=future_days)
    forecast = model.predict(future)
    return model, forecast

def plot_sales_forecast(model, forecast):
    fig = model.plot(forecast)
    plt.title("📈 Retail Sales Forecast")
    return fig
