import React, { useEffect, useState, useRef } from 'react';
import { useParams } from 'react-router-dom';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { useEnvironment } from "../context/EnvironmentContext";
import { useModel } from "../context/ModelContext";
import ChatInputBar from './shared/ChatInputBar';
import { io, Socket } from 'socket.io-client';
import LoadingSpinner from './shared/LoadingSpinner';

let mountCount = 0;

const ChatSession = () => {
  const { sessionId } = useParams<{ sessionId: string }>();
  const [messages, setMessages] = useState<any[]>([]);
  const [loading, setLoading] = useState(true);
  const [sending, setSending] = useState(false);
  const [input, setInput] = useState('');
  const [error, setError] = useState<string | null>(null);
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const socketRef = useRef<Socket | null>(null);
  const { selectedEnv } = useEnvironment();
  const { selectedModel } = useModel();

  console.log("ChatSession mounted", sessionId);

  useEffect(() => {
    mountCount += 1;
    console.log("ChatSession useEffect mount count:", mountCount);
  }, []);

  useEffect(() => {
    // Only create the socket once
    if (!socketRef.current) {
      socketRef.current = io('/', { withCredentials: true });
    }
    return () => {
      socketRef.current?.disconnect();
    };
  }, []);

  useEffect(() => {
    if (!socketRef.current || !sessionId) return;

    // Join the room
    console.log('[SocketIO] Emitting join_chat_session', sessionId);
    socketRef.current.emit('join_chat_session', { session_id: sessionId });

    // Handler for new_message
    const handler = (msg: any) => {
      console.log('[SocketIO] Received new_message:', msg);
      setMessages(prev => [...prev, msg]);
    };
    socketRef.current.on('new_message', handler);

    // Cleanup handler only (not the socket)
    return () => {
      socketRef.current?.off('new_message', handler);
    };
  }, [sessionId]);

  useEffect(() => {
    if (!sessionId) return;
    setLoading(true);
    fetch(`/api/chat_sessions/${sessionId}/messages`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        setMessages(data.messages || []);
      })
      .catch(err => {
        setMessages([]);
        setError("Failed to load messages.");
        console.error("Error loading messages:", err);
      })
      .finally(() => setLoading(false));
  }, [sessionId]);

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
    setError(null);
    // Map selectedModel to workflow_id
    // 1: ToxIndex RAP probra_task
    // 2: ToxIndex Vanilla plain_openai_task
    // 3: ToxIndex Pathway openai_json_schema_task
    // 4: ToxIndex 4th probra_task
    // 5: ToxIndex 5th probra_task
    let workflow_id = 1;
    if (selectedModel === "toxindex-rap") workflow_id = 1;
    else if (selectedModel === "toxindex-vanilla") workflow_id = 2;
    else if (selectedModel === "toxindex-pathway") workflow_id = 3;
    else if (selectedModel === "toxindex-4th" || selectedModel === "toxindex-5th") {
      setError("This workflow is not yet supported.");
      setSending(false);
      return;
    }
    try {
      const res = await fetch('/api/tasks', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          message: input.trim(),
          workflow: workflow_id,
          environment_id: selectedEnv || undefined,
          sid: sessionId,
        }),
        credentials: 'include'
      });
      const data = await res.json();
      if (data.error) {
        setError(data.error);
      } else {
        setMessages(prev => [...prev, { content: input.trim(), role: 'user' }]);
        setInput('');
      }
    } catch (error) {
      setError('Error sending message.');
      console.error('Error sending message:', error);
    } finally {
      setSending(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col flex-1" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)' }}>
      <div className="flex-1 flex flex-col max-w-5xl mx-auto w-full px-4 pt-20 py-0">
        {/* Chat history box */}
        <div
          className="flex-1 overflow-y-auto w-full max-w-5xl mx-auto px-20 py-0 bg-transparent rounded-lg shadow-inner"
          style={{
            scrollBehavior: 'smooth',
            background: 'transparent',
            minHeight: 0, // allow flex child to shrink
            maxHeight: '75vh', // limit height so only this box scrolls
          }}
        >
          {loading ? (
            <LoadingSpinner text="Loading messages..." />
          ) : (
            <>
              {messages.length === 0 && (
                <div className="text-gray-400 text-center mt-12">No messages yet. Start the conversation!</div>
              )}
              {messages.map((msg, idx) => (
                <div key={idx} className={`flex mb-4 ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`} style={{ width: '100%' }}>
                  {msg.role === 'user' ? (
                    <div
                      className="rounded-2xl px-5 py-3 shadow-md max-w-[70%] bg-gradient-to-br from-purple-600 to-green-500 text-white ml-auto"
                      style={{
                        textAlign: 'right',
                        fontSize: '1rem',
                        wordBreak: 'break-word',
                        borderBottomRightRadius: 8,
                        borderBottomLeftRadius: 24,
                        borderTopLeftRadius: 24,
                        borderTopRightRadius: 24,
                        marginLeft: 32,
                        marginRight: 0,
                        boxShadow: '0 2px 8px 0 rgba(34,197,94,0.10)'
                      }}
                    >
                      <ReactMarkdown remarkPlugins={[remarkGfm]}>{msg.content}</ReactMarkdown>
                    </div>
                  ) : (
                    <div
                      className="text-gray-100 w-full"
                      style={{
                        textAlign: 'left',
                        fontSize: '1rem',
                        wordBreak: 'break-word',
                        marginLeft: 0,
                        lineHeight: 1.8,
                      }}
                    >
                      <ReactMarkdown
                        remarkPlugins={[remarkGfm]}
                        components={{
                          td: ({node, ...props}) => <td className="markdown-table-cell" {...props} />
                        }}
                      >
                        {msg.content}
                      </ReactMarkdown>
                    </div>
                  )}
                </div>
              ))}
              {messages.length > 0 && messages[messages.length-1]?.role === 'user' && (
                <div className="flex justify-start mb-2">
                  <div className="rounded-xl px-4 py-2 bg-gray-700 text-white max-w-[70%] text-sm opacity-80">
                    Task is running...
                  </div>
                </div>
              )}
              <div ref={messagesEndRef} />
            </>
          )}
        </div>

        <ChatInputBar
          value={input}
          onChange={setInput}
          onSubmit={handleSubmit}
          uploading={sending}
          error={error}
        />
      </div>
    </div>
  );
};

export default ChatSession; 