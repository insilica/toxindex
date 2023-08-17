from flask import Flask
from flask_wtf import FlaskForm
from wtforms import StringField,PasswordField
from wtforms import validators

class RegistrationForm(FlaskForm):
  email    = StringField('my email', [validators.Length(min=6, max=35)])
  password = PasswordField('password', [validators.InputRequired()])

class UpdatePasswordForm(FlaskForm):
  password = PasswordField('password', [validators.InputRequired()])
  
class ForgotPasswordForm(FlaskForm):
  email    = StringField('my email', [validators.Length(min=6, max=35)])