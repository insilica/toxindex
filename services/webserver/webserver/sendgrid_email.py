import sendgrid
import os
from sendgrid.helpers.mail import Mail, Email, To, Content
from webserver.resources import resources as Resources

class BBSendgrid():
  
  @staticmethod
  def send(email,subject,content):
    sg = sendgrid.SendGridAPIClient(api_key=Resources.get_env('SENDGRID_API_KEY'))
    
    from_email = Email("info@insilica.co")  # Change to your verified sender
    to_email = To(email)  # Change to your recipient
    
    content = Content("text/html", content)
    mail = Mail(from_email, to_email, subject, content)

    # Get a JSON-ready representation of the Mail object
    mail_json = mail.get()

    # Send an HTTP POST request to /mail/send
    response = sg.client.mail.send.post(request_body=mail_json)
    print(response.status_code)
    print(response.headers)
    return response