"""
Retail Sales Forecasting 
Author: Phat Nguyen Cong - https://github.com/paht2005

"""

import pandas as pd
import matplotlib.pyplot as plt
from prophet import Prophet
import logging
import argparse
import sys
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def generate_synthetic_sales_data(start_date="2024-01-01", periods=120):
    """Simulate sales data with weekly seasonality and some noise."""
    try:
        dates = pd.date_range(start=start_date, periods=periods)
        sales = [
            200 + (i % 7) * 10 + (i % 15) * (-1) ** i + (5 if i % 20 == 0 else 0)
            for i in range(periods)
        ]
        df = pd.DataFrame({'ds': dates, 'y': sales})
        return df
    except Exception as e:
        logging.error(f"Failed to generate synthetic data: {e}")
        sys.exit(1)

def build_and_train_model(df, daily=True):
    """Initialize and train the Prophet model."""
    try:
        model = Prophet(daily_seasonality=daily, yearly_seasonality=False)
        model.fit(df)
        return model
    except Exception as e:
        logging.error(f"Model training failed: {e}")
        sys.exit(1)

def extend_forecast(model, future_days=30):
    """Forecast future sales."""
    try:
        future_df = model.make_future_dataframe(periods=future_days)
        forecast = model.predict(future_df)
        return forecast
    except Exception as e:
        logging.error(f"Forecast generation failed: {e}")
        sys.exit(1)

def plot_forecast(model, forecast, output_path=None):
    """Plot the forecast and optionally save the plot."""
    try:
        fig = model.plot(forecast)
        plt.title("Retail Sales Forecast")
        if output_path:
            plt.savefig(output_path, bbox_inches='tight')
            logging.info(f"Plot saved to {output_path}")
        else:
            plt.show()
    except Exception as e:
        logging.error(f"Failed to plot forecast: {e}")
        sys.exit(1)

def main(args):
    df = generate_synthetic_sales_data()
    model = build_and_train_model(df)
    forecast = extend_forecast(model, future_days=args.days)
    plot_forecast(model, forecast, output_path=args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retail Sales Forecasting Tool")
    parser.add_argument('--days', type=int, default=30, help='Number of days to forecast')
    parser.add_argument('--output', type=str, help='Optional output path for saving the plot')
    args = parser.parse_args()

    main(args)
