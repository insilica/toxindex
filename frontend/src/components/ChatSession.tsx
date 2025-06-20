import React, { useEffect, useState, useRef } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { FaArrowLeft } from 'react-icons/fa';
import EnvironmentSelector from './shared/EnvironmentSelector';
import LoadingSpinner from './shared/LoadingSpinner';
import { useEnvironmentSelection } from './shared/utils/useEnvironmentSelection';

interface Environment {
  environment_id: string;
  title: string;
}

interface ChatSessionProps {
  environments: Environment[];
  selectedEnv?: string;
  setSelectedEnv?: (envId: string) => void;
  loadingEnvironments?: boolean;
}

const ChatSession: React.FC<ChatSessionProps> = ({
  environments,
  selectedEnv,
  setSelectedEnv,
  loadingEnvironments = false
}) => {
  const { session_id } = useParams<{ session_id: string }>();
  const navigate = useNavigate();
  const [messages, setMessages] = useState<any[]>([]);
  const [loading, setLoading] = useState(true);
  const [sending, setSending] = useState(false);
  const [input, setInput] = useState('');
  const messagesEndRef = useRef<HTMLDivElement>(null);

  useEnvironmentSelection(environments, selectedEnv, setSelectedEnv);

  useEffect(() => {
    if (!session_id) return;
    setLoading(true);
    fetch(`/api/chat-sessions/${session_id}`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        setMessages(data.messages || []);
        setLoading(false);
      });
  }, [session_id]);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!input.trim() || sending) return;

    setSending(true);
    try {
      const res = await fetch(`/api/chat-sessions/${session_id}/messages`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ content: input.trim() }),
        credentials: 'include'
      });
      const data = await res.json();
      if (data.success) {
        setMessages(prev => [...prev, data.message]);
        setInput('');
      }
    } catch (error) {
      console.error('Error sending message:', error);
    } finally {
      setSending(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col flex-1" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)' }}>
      <div className="flex-1 flex flex-col max-w-4xl mx-auto w-full px-4 py-8">
        <div className="flex items-center mb-8">
          <button
            onClick={() => navigate(-1)}
            className="mr-4 text-gray-400 hover:text-white transition-colors"
          >
            <FaArrowLeft size={24} />
          </button>
          <EnvironmentSelector
            environments={environments}
            selectedEnv={selectedEnv}
            onEnvironmentChange={setSelectedEnv}
            loadingEnvironments={loadingEnvironments}
            variant="large"
          />
        </div>

        <div className="flex-1 bg-black bg-opacity-40 rounded-lg p-6 mb-8 overflow-y-auto">
          {loading ? (
            <div className="flex-1 flex items-center justify-center">
              <LoadingSpinner size="large" />
            </div>
          ) : messages.length === 0 ? (
            <div className="text-gray-400 text-center">No messages yet. Start the conversation!</div>
          ) : (
            <div className="space-y-6">
              {messages.map((msg, idx) => (
                <div key={idx} className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}>
                  <div className={`max-w-3xl p-4 rounded-lg ${msg.role === 'user' ? 'bg-green-900 bg-opacity-40' : 'bg-gray-800 bg-opacity-40'}`}>
                    <pre className="whitespace-pre-wrap font-sans text-white">{msg.content}</pre>
                  </div>
                </div>
              ))}
              <div ref={messagesEndRef} />
            </div>
          )}
        </div>

        <form onSubmit={handleSubmit} className="flex gap-4">
          <input
            type="text"
            value={input}
            onChange={e => setInput(e.target.value)}
            placeholder="Type your message..."
            className="flex-1 bg-black bg-opacity-40 text-white px-4 py-3 rounded-lg focus:outline-none focus:ring-2 focus:ring-green-500"
            disabled={sending}
          />
          <button
            type="submit"
            className={`px-6 py-3 bg-green-600 text-white rounded-lg font-medium hover:bg-green-700 transition-colors ${sending ? 'opacity-50 cursor-not-allowed' : ''}`}
            disabled={sending}
          >
            {sending ? <LoadingSpinner size="small" text="" /> : 'Send'}
          </button>
        </form>
      </div>
    </div>
  );
};

export default ChatSession; 