import stripe
import os
import ssl
import requests
from requests.adapters import HTTPAdapter

# Configure Stripe
stripe.api_version = '2020-08-27'
stripe.api_key = os.getenv('STRIPE_SECRET_KEY')

# Create SSL context with TLS 1.2 requirement
context = ssl.create_default_context()
context.minimum_version = ssl.TLSVersion.TLSv1_2

# Create a session with our SSL context
session = requests.Session()
adapter = HTTPAdapter()
adapter.init_poolmanager(connections=10, maxsize=10, ssl_context=context)
session.mount('https://', adapter)

# Configure Stripe to use our session
stripe.default_http_client = stripe.http_client.RequestsClient(verify_ssl_certs=True, session=session)

def create_customer(email):
  return stripe.Customer.create(email=email)

def delete_customer(stripe_customer_id):
  return stripe.Customer.delete(stripe_customer_id)

def create_customer_portal_session(stripe_customer_id):
  return stripe.billing_portal.Session.create(customer=stripe_customer_id)
