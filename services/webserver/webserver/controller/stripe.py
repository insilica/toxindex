import stripe
import os

class ToxindexStripe():
    
  stripe.api_version = '2020-08-27'
  stripe.api_key = os.getenv('STRIPE_TEST_SECRET_KEY')

  @staticmethod
  def create_customer(email):
    return stripe.Customer.create(email=email)

  @staticmethod
  def delete_customer(stripe_customer_id):
    return stripe.Customer.delete(stripe_customer_id)

  @staticmethod
  def create_customer_portal_session(stripe_customer_id):
    return stripe.billing_portal.Session.create(customer=stripe_customer_id)
