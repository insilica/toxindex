import os
import json
import openai

openai.api_key = os.environ.get('OPENAI_API_KEY')

def generate_response(message, conversation_history=None):
    """
    Generate an AI response using an external API service
    
    Args:
        message (str): The user's message
        conversation_history (list): Optional list of previous messages
        
    Returns:
        str: The AI response
    """
    # For a real integration, you'd use something like this:
    # 
    # api_key = os.environ.get('OPENAI_API_KEY')
    # url = "https://api.openai.com/v1/chat/completions"
    # 
    # headers = {
    #     "Content-Type": "application/json",
    #     "Authorization": f"Bearer {api_key}"
    # }
    # 
    # # Format conversation history
    # messages = []
    # if conversation_history:
    #     for msg in conversation_history:
    #         messages.append({"role": msg["role"], "content": msg["content"]})
    # 
    # # Add the current message
    # messages.append({"role": "user", "content": message})
    # 
    # data = {
    #     "model": "gpt-3.5-turbo",
    #     "messages": messages,
    #     "temperature": 0.7
    # }
    # 
    # response = requests.post(url, headers=headers, data=json.dumps(data))
    # response_data = response.json()
    # 
    # return response_data["choices"][0]["message"]["content"]
    
    # For demo purposes, return a simple response
    return f"You said: {message}. This is a placeholder response. In a real application, this would come from an AI model."


def generate_title(message: str) -> str:
    """Generate a short title summarizing the given message using OpenAI."""
    try:
        resp = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "Provide a short title summarizing the user's request."},
                {"role": "user", "content": message},
            ],
            max_tokens=10,
            temperature=0.5,
        )
        return resp.choices[0].message["content"].strip()
    except Exception:
        # Fallback title
        return message.strip()[:20]
