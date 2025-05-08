import stripe
import os
  
stripe.api_version = '2020-08-27'
stripe.api_key = os.getenv('STRIPE_TEST_SECRET_KEY')

def create_customer(email):
  return stripe.Customer.create(email=email)

def delete_customer(stripe_customer_id):
  return stripe.Customer.delete(stripe_customer_id)

def create_customer_portal_session(stripe_customer_id):
  return stripe.billing_portal.Session.create(customer=stripe_customer_id)
