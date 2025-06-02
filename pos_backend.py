# POS Interface (Backend Only)
# Phat Nguyen Cong - https://github.com/paht2005

import pandas as pd
import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple

# Setup basic logging
logging.basicConfig(level=logging.INFO)

@dataclass
class Product:
    id: str
    name: str
    price: float

# Simulated Product DB
PRODUCTS_DB: Dict[str, Product] = {
    "001": Product(id="001", name="Paracetamol", price=1.5),
    "002": Product(id="002", name="Ibuprofen", price=2.0),
    "003": Product(id="003", name="ORS Pack", price=0.8)
}

def fetch_product_by_id(pid: str) -> Product:
    """Fetch product from the database by its ID."""
    product = PRODUCTS_DB.get(pid)
    if not product:
        logging.warning(f"Product ID '{pid}' not found in the database.")
    return product

def calculate_subtotal(price: float, quantity: int) -> float:
    """Calculate subtotal for a single product."""
    return round(price * quantity, 2)

def create_invoice(cart_items: Dict[str, int]) -> Tuple[List[Dict], float]:
    """
    Generate an invoice from the shopping cart.

    Args:
        cart_items: Dictionary mapping product IDs to quantities.

    Returns:
        Tuple of (line items list, total price).
    """
    invoice_lines = []
    total_price = 0.0

    for pid, qty in cart_items.items():
        product = fetch_product_by_id(pid)
        if product:
            line_subtotal = calculate_subtotal(product.price, qty)
            line_item = {
                "product": product.name,
                "qty": qty,
                "unit": product.price,
                "subtotal": line_subtotal
            }
            invoice_lines.append(line_item)
            total_price += line_subtotal
        else:
            logging.info(f"Skipping unknown product ID: {pid}")

    total_price = round(total_price, 2)
    return invoice_lines, total_price
